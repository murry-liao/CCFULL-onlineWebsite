#!/usr/bin/env python3
"""Discovery and parsing of the local Expdata fusion cross-section library."""

from __future__ import annotations

from pathlib import Path
import re


from physical_model.paths import EXPERIMENTAL_DATA_ROOT


EXPDATA_ROOT = EXPERIMENTAL_DATA_ROOT
PATTERN = re.compile(
    r"^fusion-exp-(?P<z1>\d+)\.(?P<a1>\d+)\.(?P<z2>\d+)\.(?P<a2>\d+)"
    r"-(?P<kind>[^-]+)-(?P<year>\d+)-(?P<dataset>\d+)\.txt$"
)
BARRIER_PATTERN = re.compile(
    r"^fusion-barrier-(?P<z1>\d+)\.(?P<a1>\d+)\.(?P<z2>\d+)\.(?P<a2>\d+)"
    r"-(?P<label>.+)\.txt$"
)


def dataset_metadata(path: str | Path) -> dict:
    path = Path(path)
    match = PATTERN.match(path.name)
    if not match:
        raise ValueError(f"unrecognized Expdata filename: {path.name}")
    values = match.groupdict()
    result = {key: int(values[key]) for key in ("z1", "a1", "z2", "a2", "year", "dataset")}
    result.update(kind=values["kind"], path=path.resolve())
    result["label"] = (
        f"Z={result['z1']}, A={result['a1']} + Z={result['z2']}, A={result['a2']} · "
        f"{result['kind']} · {result['year']} · dataset {result['dataset']}"
    )
    return result


def discover_datasets(apro: int, zpro: int, atar: int, ztar: int) -> list[dict]:
    requested = ((int(zpro), int(apro)), (int(ztar), int(atar)))
    found = []
    if not EXPDATA_ROOT.exists():
        return found
    for path in EXPDATA_ROOT.rglob("fusion-exp-*.txt"):
        try:
            metadata = dataset_metadata(path)
        except ValueError:
            continue
        stored = ((metadata["z1"], metadata["a1"]), (metadata["z2"], metadata["a2"]))
        if stored == requested or stored == requested[::-1]:
            found.append(metadata)
    return sorted(found, key=lambda item: (item["year"], item["kind"], item["dataset"]))


def load_dataset(path: str | Path) -> list[dict]:
    metadata = dataset_metadata(path)
    rows = []
    for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 5:
            raise ValueError(f"{path}:{line_number}: expected five columns")
        energy_cm, cross_section, uncertainty_lower, uncertainty_upper, energy_lab = map(float, fields)
        rows.append({
            "energy_cm_mev": energy_cm,
            "fusion_cross_section_mb": cross_section,
            "uncertainty_lower_mb": uncertainty_lower,
            "uncertainty_upper_mb": uncertainty_upper,
            "energy_lab_mev": energy_lab,
            "dataset": metadata["label"],
            "source_file": str(metadata["path"]),
        })
    if not rows:
        raise ValueError(f"experimental dataset is empty: {path}")
    return rows


def discover_barrier_datasets(apro: int, zpro: int, atar: int, ztar: int) -> list[dict]:
    requested = ((int(zpro), int(apro)), (int(ztar), int(atar)))
    found = []
    if not EXPDATA_ROOT.exists():
        return found
    for path in EXPDATA_ROOT.rglob("fusion-barrier-*.txt"):
        match = BARRIER_PATTERN.match(path.name)
        if not match:
            continue
        values = match.groupdict()
        stored = ((int(values["z1"]), int(values["a1"])),
                  (int(values["z2"]), int(values["a2"])))
        if stored == requested or stored == requested[::-1]:
            found.append({
                "path": path.resolve(),
                "label": values["label"].replace("-", " "),
            })
    return sorted(found, key=lambda item: item["label"])


def load_barrier_dataset(path: str | Path) -> list[dict]:
    path = Path(path)
    match = BARRIER_PATTERN.match(path.name)
    if not match:
        raise ValueError(f"unrecognized fusion-barrier filename: {path.name}")
    rows = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"{path}:{line_number}: expected four columns")
        energy_lab, energy_cm, distribution, uncertainty = map(float, fields)
        rows.append({
            "energy_cm_mev": energy_cm,
            "energy_lab_mev": energy_lab,
            "barrier_distribution_mb_per_mev": distribution,
            "uncertainty_mb_per_mev": uncertainty,
            "dataset": match.group("label").replace("-", " "),
            "source_file": str(path.resolve()),
        })
    if not rows:
        raise ValueError(f"fusion-barrier dataset is empty: {path}")
    return rows
