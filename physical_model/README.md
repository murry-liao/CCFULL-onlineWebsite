# Physical model

- `solvers/standard_numerov/ccfull.f90` is the high-fidelity model run by the
  website.
- `solvers/hagino_original/ccfull.for` is the retained original CCFULL baseline.
- `interface/physical_runner.py` builds both programs, runs them in a temporary
  directory, and parses their cross sections.
- `interface/fortran_compiler.py` selects the gfortran that generates code for
  the running Python's architecture.
- `../data/input_templates/` contains physical CCFULL input records.
- `../data/experimental/4-final-ngt4/` is the active experimental library: 410
  curated fusion datasets covering 369 reaction systems.
- `build/` is generated and may be recreated from source.
