# Experimental data and input templates

This directory contains inputs of the website. It contains no executable code.

- `experimental/4-final-ngt4/`: curated fusion excitation-function data.
- `input_templates/`: editable CCFULL input templates.

The active experimental library contains 410 excitation-function files for
369 reaction systems and 7,311 measured points. Each row stores
`E_cm, sigma_fusion, uncertainty_lower, uncertainty_upper, E_lab`.

`input_templates/o16_sm154.inp` is the physical-model template. The website
loads it to fill in the default parameters, so editing this file changes the
defaults shown on the page.
