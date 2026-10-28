

# CCFULL Online Calculation Website

A Streamlit-based web platform for online CCFULL fusion reaction calculations. The site runs the **standard Numerov solver (Fortran 90)** alongside **Professor Hagino's original CCFULL (Fortran 77)** as a baseline, allows direct comparison of results with experimental data, and exports the barrier distribution.

---

## ✨ Features

- 🔬 Dual solvers: standard Numerov solver + original CCFULL baseline
- 📊 Direct comparison of results with experimental data
- 📈 Barrier distribution export
- 🧩 Curated experimental datasets (from the NRV website) and CCFULL input templates
- ⚡ Automatic compilation of Fortran programs on first run, with binaries reused thereafter

---

## 📁 Project Structure

```text
CCFULL_Website/
├── app.py                Streamlit page
├── start.py              launcher
├── requirements.txt      dependency list
├── physical_model/
│   ├── solvers/          standard Numerov and Hagino original source
│   ├── interface/        input, compiler, solver-runner, and data adapters
│   └── build/            generated binaries
└── data/
    ├── experimental/     curated 4-final-ngt4 fusion datasets from NRV websites
    └── input_templates/  CCFULL inputs offered by the page
```

---

## 🛠️ Requirements

- Python **3.11** or newer
- **gfortran** with LAPACK support
- Python packages listed in `requirements.txt`

---

## 🚀 Installation & Startup

### 1. Prepare the environment

```bash
conda activate ccfull-rom
```

> For a fresh environment, install the dependencies first:

```bash
python -m pip install -r requirements.txt
```

### 2. Start the website

```bash
python start.py
```

Then open <http://127.0.0.1:8501/> in your browser.

---

## 📝 Usage Notes

- **First calculation**: Both Fortran programs are compiled into `physical_model/build/` using the gfortran that targets the running Python's architecture.
- **Subsequent calculations**: The compiled binaries are reused automatically — no recompilation needed.
- **About the baseline curve**: Hagino's original source reads the legacy four-field rotor record, so its curve corresponds to the `beta6 = 0` baseline whenever `Beta6T` is non-zero.

---

## 📄 License

<!-- Add license information here, e.g. MIT / Apache-2.0 -->

---

## 🙏 Acknowledgements

- Original CCFULL program: Prof. K. Hagino
- Experimental datasets: NRV  website

---

> 
