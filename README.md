

# CCFULL Online Calculation Website

A Streamlit-based web platform for online CCFULL fusion reaction calculations. The site runs the **standard Numerov solver (Fortran 90)** alongside **Professor Hagino's original CCFULL (Fortran 77)** as a baseline, allows direct comparison of results with experimental data, and plots the barrier distribution.

---

## ✨ Features

- 🔬 Dual solvers: standard Numerov solver + original CCFULL baseline
- 📊 Direct comparison of results with experimental fusion cross sections
- 📈 Barrier distribution D(E) of the calculation and of the experimental data
- 💾 CSV download of the cross sections and of the CCFULL input file that was used
- 🧩 Curated experimental datasets (from the NRV website) and a default CCFULL input template
- ⚡ Automatic compilation of Fortran programs on first run, with binaries reused thereafter

---

## 📁 Project Structure

```text
CCFULL-onlineWebsite/
├── app.py                Streamlit page
├── start.py              launcher
├── requirements.txt      Python dependencies
├── physical_model/
│   ├── solvers/          standard Numerov and Hagino original source
│   ├── interface/        input, compiler, solver-runner, and data adapters
│   └── build/            generated binaries (created on first run)
└── data/
    ├── experimental/     curated 4-final-ngt4 fusion datasets from NRV websites
    └── input_templates/  default CCFULL input (16O + 154Sm) loaded by the page
```

---

## 🛠️ Requirements

- **macOS or Linux** — on Windows use WSL (the solver runner relies on `fcntl`, which native Windows lacks)
- Python **3.11** or newer
- **gfortran** with LAPACK and BLAS (the Fortran programs are linked with `-llapack -lblas`)
- Python packages listed in `requirements.txt` (`streamlit`, `pandas`, `numpy`, `plotly`)

System packages:

| Platform | Command |
| --- | --- |
| macOS (Homebrew) | `brew install gcc` — provides `gfortran`; LAPACK/BLAS come from Apple's Accelerate framework |
| Debian / Ubuntu | `sudo apt install gfortran liblapack-dev libblas-dev` |

Developed and tested on macOS (Apple silicon) with Python 3.12, Streamlit 1.66, pandas 3.0–3.1, numpy 2.5, plotly 7.1 and GNU Fortran 15.1 (Homebrew GCC).

---

## 🚀 Installation & Startup

### 1. Prepare the environment

Create a Python environment and install the dependencies (conda is shown; `python -m venv .venv` works as well):

```bash
conda create -n ccfull-website python=3.12
conda activate ccfull-website
python -m pip install -r requirements.txt
```

Check that a native `gfortran` is on your `PATH`:

```bash
gfortran --version
```

### 2. Start the website

```bash
python start.py
```

Equivalently, run `streamlit run app.py`. Then open <http://127.0.0.1:8501/> in your browser (if port 8501 is busy, Streamlit uses the next free port — use the address printed in the terminal).

---

## 📖 Using the Website

1. Fill in the **Parameter settings** panel: reaction system and radii, coupling parameters, then the nuclear potential and numerical grid. The defaults are loaded from `data/input_templates/o16_sm154.inp` (16O + 154Sm).
2. Optionally tick **Include experimental fusion cross-section data** to overlay the datasets that match the projectile and target, and set the energy step dE used for the experimental barrier distribution.
3. Click **Run Physical Model**. The Numerov solver and Hagino's original CCFULL are both run.
4. Read the cross-section and D(E) plots, the timing and baseline-comparison metrics and the results table; download the CSV files or the `ccfull.inp` that was used.

---

## 📝 Usage Notes

- **First calculation**: Both Fortran programs are compiled into `physical_model/build/` using the gfortran that targets the running Python's architecture.
- **Subsequent calculations**: The compiled binaries are reused automatically — no recompilation needed.
- **About the baseline curve**: Hagino's original source reads the legacy four-field rotor record, so its curve corresponds to the `beta6 = 0` baseline whenever `Beta6T` is non-zero.
- **Troubleshooting**:
  - `no gfortran on PATH generates arm64 code for this Python` — on Apple silicon an Intel `gfortran` (for example in `/usr/local/bin`) can shadow the native one. Install a native one (`brew install gcc`) or put it first on `PATH`.
  - Linker errors mentioning `-llapack` or `-lblas` — install the LAPACK/BLAS development packages (see Requirements).
  - The Fortran programs are rebuilt automatically when a solver source changes; to force a rebuild, delete `physical_model/build/`.

---

## 📄 License

<!-- Add license information here, e.g. MIT / Apache-2.0 -->

---

## 🙏 Acknowledgements

- Original CCFULL program: Prof. K. Hagino — K. Hagino, N. Rowley and A. T. Kruppa, *Comput. Phys. Commun.* **123**, 143 (1999)
- Experimental datasets: NRV website
