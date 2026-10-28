#!/usr/bin/env python3
"""Generic input and execution interface for the real-potential CCFULL model."""

from __future__ import annotations

import math
from pathlib import Path
import re


def _items(line: str) -> list[str]:
    return [item for item in re.split(r"[,\s]+", line.strip()) if item]


def clean_lines(path: str | Path) -> list[str]:
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if len(lines) < 9:
        raise ValueError("CCFULL input must contain nine data records")
    return lines[:9]


def parse_input(path: str | Path) -> dict:
    lines = clean_lines(path)
    reaction = list(map(float, _items(lines[0])))
    geometry = list(map(float, _items(lines[1])))
    target1 = list(map(float, _items(lines[2])))
    target2 = list(map(float, _items(lines[3])))
    projectile = list(map(float, _items(lines[4])))
    transfer = list(map(float, _items(lines[5])))
    potential = list(map(float, _items(lines[6])))
    energy = list(map(float, _items(lines[7])))
    radial = list(map(float, _items(lines[8])))
    if any(len(row) < size for row, size in (
        (reaction, 4), (geometry, 4), (target1, 4), (target2, 4),
        (projectile, 4), (transfer, 2), (potential, 3), (energy, 3), (radial, 2)
    )):
        raise ValueError("one or more CCFULL input records have too few values")
    target_mode = int(geometry[3])
    if target_mode == 1:
        # Extended rotor record: E2T,beta2,beta4,beta6,NrotT.  Preserve
        # compatibility with the original four-field record by inserting a
        # zero static beta6 before NrotT.
        target1 = (target1[:5] if len(target1) >= 5 else
                   [target1[0], target1[1], target1[2], 0.0, target1[3]])
    else:
        target1 = target1[:4]
    return {
        "apro": int(reaction[0]), "zpro": int(reaction[1]),
        "atar": int(reaction[2]), "ztar": int(reaction[3]),
        "r0p": geometry[0], "projectile_mode": int(geometry[1]),
        "r0t": geometry[2], "target_mode": target_mode,
        "target1": target1, "target2": target2[:4],
        "projectile": projectile[:4],
        "qtrans": transfer[0], "ftr": transfer[1],
        "v0": potential[0], "r0": potential[1], "a0": potential[2],
        "emin": energy[0], "emax": energy[1], "de": energy[2],
        "rmax": radial[0], "dr": radial[1],
        "lmax": 50, "tail_window": 2, "tail_tolerance": 1e-4,
    }


def _line(values) -> str:
    return ",".join(str(int(value)) if isinstance(value, int) else f"{float(value):.12g}"
                    for value in values)


def validate_config(config: dict) -> None:
    for key in ("apro", "zpro", "atar", "ztar"):
        if int(config[key]) <= 0:
            raise ValueError(f"{key} must be positive")
    if int(config["zpro"]) > int(config["apro"]) or int(config["ztar"]) > int(config["atar"]):
        raise ValueError("nuclear charge cannot exceed mass number")
    if int(config["projectile_mode"]) not in (-1, 0, 1) or int(config["target_mode"]) not in (-1, 0, 1):
        raise ValueError("coupling mode must be -1 (inert), 0 (vibration), or 1 (rotation)")
    for key in ("r0p", "r0t", "v0", "r0", "a0", "de", "rmax", "dr"):
        if float(config[key]) <= 0:
            raise ValueError(f"{key} must be positive")
    if float(config["emax"]) < float(config["emin"]):
        raise ValueError("emax must be at least emin")
    if not 0 <= int(config["lmax"]) <= 200:
        raise ValueError("lmax must be in [0, 200]")
    if int(config["tail_window"]) < 1 or float(config["tail_tolerance"]) <= 0:
        raise ValueError("tail controls must be positive")
    if energy_values(config)[-1] > float(config["emax"]) + 1e-8:
        raise ValueError("invalid energy grid")


def input_text(config: dict) -> str:
    validate_config(config)
    rows = [
        _line([int(config["apro"]), int(config["zpro"]), int(config["atar"]), int(config["ztar"])]),
        _line([config["r0p"], int(config["projectile_mode"]), config["r0t"], int(config["target_mode"])]),
        _line(config["target1"]), _line(config["target2"]), _line(config["projectile"]),
        _line([config["qtrans"], config["ftr"]]),
        _line([config["v0"], config["r0"], config["a0"]]),
        _line([config["emin"], config["emax"], config["de"]]),
        _line([config["rmax"], config["dr"]]),
    ]
    return "\n".join(rows) + "\n"


def energy_values(config: dict) -> list[float]:
    first, last, step = map(float, (config["emin"], config["emax"], config["de"]))
    if step <= 0 or last < first:
        raise ValueError("energy grid requires de > 0 and emax >= emin")
    count = int(math.floor((last - first) / step + 1e-9))
    values = [first + index * step for index in range(count + 1)]
    if len(values) > 401:
        raise ValueError("energy grid is limited to 401 points in the web interface")
    return values
