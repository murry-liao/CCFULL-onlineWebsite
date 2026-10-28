"""Stable paths for the physical-model package."""

from pathlib import Path


PHYSICAL_MODEL_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PHYSICAL_MODEL_ROOT.parent
DATA_ROOT = PROJECT_ROOT / "data"
INPUT_TEMPLATE_ROOT = DATA_ROOT / "input_templates"
STANDARD_NUMEROV_SOURCE = (
    PHYSICAL_MODEL_ROOT / "solvers" / "standard_numerov" / "ccfull.f90")
HAGINO_SOURCE = (
    PHYSICAL_MODEL_ROOT / "solvers" / "hagino_original" / "ccfull.for")
DEFAULT_INPUT = INPUT_TEMPLATE_ROOT / "o16_sm154.inp"
EXPERIMENTAL_DATA_ROOT = (
    DATA_ROOT / "experimental" / "4-final-ngt4")
BUILD_ROOT = PHYSICAL_MODEL_ROOT / "build"
