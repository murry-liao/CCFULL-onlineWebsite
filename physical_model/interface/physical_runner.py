#!/usr/bin/env python3
"""Build and run selectable high-fidelity CCFULL implementations."""

from __future__ import annotations

import csv
import fcntl
from pathlib import Path
import subprocess
import tempfile
import time

from .fortran_compiler import gfortran, machine
from .input_parser import input_text
from physical_model.paths import BUILD_ROOT, HAGINO_SOURCE, STANDARD_NUMEROV_SOURCE


BUILD = BUILD_ROOT / "standalone"
CODE_SPECS = {
    "real_numerov": {
        "label": "Real potential · standard Numerov",
        "source": STANDARD_NUMEROV_SOURCE,
        "method": "original_program",
        "observable": "Fusion cross section (IWBC)",
    },
    "hagino_original": {
        "label": "Hagino original CCFULL · Fortran 77 baseline",
        "source": HAGINO_SOURCE,
        "method": "hagino_original",
        "legacy": True,
        "observable": "Fusion cross section (IWBC)",
    },
}


def build_original(code: str) -> Path:
    spec = CODE_SPECS[code]
    compiler = gfortran()
    BUILD.mkdir(parents=True, exist_ok=True)
    binary = BUILD / f"{code}_{machine()}"
    source = Path(spec["source"])
    newest_source = max(source.stat().st_mtime, Path(__file__).stat().st_mtime)
    if binary.exists() and binary.stat().st_mtime >= newest_source:
        return binary
    with (BUILD / f".{code}.build.lock").open("w") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if binary.exists() and binary.stat().st_mtime >= newest_source:
            return binary
        command = [compiler, "-O2"]
        if spec.get("legacy"):
            command.extend(("-std=legacy", "-ffixed-line-length-none"))
        else:
            command.append("-ffree-line-length-none")
        command.extend(("-o", str(binary), str(source)))
        if not spec.get("legacy"):
            command.extend(("-llapack", "-lblas"))
        completed = subprocess.run(
            command,
            cwd=BUILD, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        if completed.returncode:
            raise RuntimeError(f"failed to build {code}:\n{completed.stdout}")
    return binary


def _parse_sigma(path: Path, cross_section_key: str) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines()[1:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        energy, cross_section, mean_l, p0 = map(float, fields[:4])
        rows.append({
            "energy_mev": energy,
            cross_section_key: cross_section,
            "mean_angular_momentum_hbar": mean_l,
            "s_wave_probability": p0,
        })
    if not rows:
        raise RuntimeError("selected physical code produced no sigma.dat rows")
    return rows


def _hagino_input_text(config: dict) -> str:
    legacy_config = dict(config)
    target1 = list(config["target1"])
    if int(config["target_mode"]) == 1 and len(target1) >= 5:
        # Hagino's original source remains unchanged and reads the legacy
        # E2T,beta2,beta4,NrotT rotor record.  It therefore provides the
        # beta6=0 baseline when the extended Numerov solver uses beta6.
        legacy_config["target1"] = target1[:3] + [target1[4]]
    lines = input_text(legacy_config).splitlines()
    ntrans = int(config.get("hagino_ntrans", 0))
    lines[5] = f"{ntrans},{float(config['qtrans']):.12g},{float(config['ftr']):.12g}"
    return "\n".join(lines) + "\n"


def _parse_hagino(work: Path) -> list[dict]:
    spin_by_energy = {}
    spin_path = work / "spin.dat"
    if spin_path.exists():
        for line in spin_path.read_text().splitlines():
            fields = line.split()
            if len(fields) >= 2:
                spin_by_energy[round(float(fields[0]), 10)] = float(fields[1])
    rows = []
    for line in (work / "cross.dat").read_text().splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        energy, cross_section = map(float, fields[:2])
        rows.append({
            "energy_mev": energy,
            "fusion_cross_section_mb": cross_section,
            "mean_angular_momentum_hbar": spin_by_energy.get(round(energy, 10)),
        })
    if not rows:
        raise RuntimeError("Hagino original CCFULL produced no cross.dat rows")
    return rows


def run_selected_code(config: dict, code: str, timeout: int = 1800):
    if code not in CODE_SPECS:
        raise ValueError(f"unknown physical code: {code}")
    spec = CODE_SPECS[code]
    binary = build_original(code).resolve()
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix=f"ccfull-{code}-") as temporary:
        work = Path(temporary)
        used_input = (_hagino_input_text(config) if spec["method"] == "hagino_original"
                      else input_text(config))
        (work / "ccfull.inp").write_text(used_input)
        if spec.get("optical"):
            optical = (float(config.get("wv0", 50.0)), float(config.get("rw0", 1.0)),
                       float(config.get("aw0", 0.3)))
            if optical[0] < 0 or optical[1] <= 0 or optical[2] <= 0:
                raise ValueError("optical parameters require Wv >= 0, rw > 0, and aw > 0")
            (work / "optical.inp").write_text(",".join(f"{value:.12g}" for value in optical) + "\n")
        completed = subprocess.run(
            [str(binary)], cwd=work,
            input="n\n" * (96 if spec["method"] == "hagino_original" else 32), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout,
        )
        expected_output = work / ("cross.dat" if spec["method"] == "hagino_original" else "sigma.dat")
        if completed.returncode or not expected_output.exists():
            raise RuntimeError(f"{spec['label']} failed:\n{completed.stdout[-8000:]}")
        if spec["method"] == "hagino_original":
            rows = _parse_hagino(work)
        else:
            cross_section_key = ("absorption_cross_section_mb" if spec.get("optical")
                                 else "fusion_cross_section_mb")
            rows = _parse_sigma(work / "sigma.dat", cross_section_key)
    return rows, [], {
        "wall_s": time.perf_counter() - started,
        "all_converged": True,
        "input": used_input,
        "code": code,
        "code_label": spec["label"],
        "source": str(spec["source"]),
        "observable": spec["observable"],
        "stdout": completed.stdout,
    }
