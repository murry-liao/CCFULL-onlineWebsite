#!/usr/bin/env python3
"""Start the CCFULL real-potential Streamlit workbench."""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
APP = ROOT / "app.py"


def main() -> int:
    try:
        import streamlit  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "Streamlit is not installed. From this directory run: "
            "python3 -m pip install -r requirements.txt"
        ) from exc
    try:
        return subprocess.call([sys.executable, "-m", "streamlit", "run", str(APP)])
    except KeyboardInterrupt:
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
