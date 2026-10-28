# CCFULL Website

Streamlit page for the online CCFULL fusion calculation. 
It runs the standard Numerov solver (fortran90) together with Professor Hagino's original CCFULL (fortran77) as a baseline.
At the same time, the result can be compared with experimental data, and exports the barrier distribution. 

```text
CCFULL_Website/
├── app.py                Streamlit page
├── start.py              launcher
├── requirements.txt
├── physical_model/
│   ├── solvers/          standard Numerov and Hagino original source
│   ├── interface/        input, compiler, solver-runner, and data adapters
│   └── build/            generated binaries
└── data/
    ├── experimental/     curated 4-final-ngt4 fusion datasets
    └── input_templates/  CCFULL inputs offered by the page
```

## Start the website

The page needs Python 3.11 or newer, a gfortran with LAPACK, and the packages
in `requirements.txt`. 

```bash
conda activate ccfull-rom
python start.py
```

In a new environment, install the packages once with
`python -m pip install -r requirements.txt`. Then open
<http://127.0.0.1:8501/>.

The first calculation compiles both Fortran programs into
`physical_model/build/` with the gfortran that generates code for the running
Python's architecture; later calculations reuse them.

Hagino's original source reads the legacy four-field rotor record, so its curve
is the `beta6 = 0` baseline whenever Beta6T is non-zero.

The `physical_model/` copy here is independent of the emulator's copy. An edit
of a solver source in one project does not reach the other.
