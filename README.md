# Herschel-Ion-Abundance-Calculator
---
## Python Scripts for Analyzing HII Sources ionic abundance.

### Overview

This set of Python files facilitates the analysis of various parameters related to HII sources, including ionic abundances, fluxes, and electron density. Below are details on how to use these scripts. To understand what is the purpose of each python file, I suggest going through the comments in the python file itself. 

### CSV File Path

The analysis relies on the following CSV file path:

- `/users/sdey/DATA/herschel.csv - Sheet2recent.csv`
- You can make your own CSV file. It should follow the structure used in the source_file.py for the pipeline to work.

### Python Files

The following Python files are utilized in the analysis:

1. `source_file.py`
2. `effecient.py`
3. `electrondensity.py`
4. `emissivity.py`
5. `freefreecalc.py`
6. `abundance_calculator.py`
7. `ionic_abundance_table.py` (output)
8. `GrayScale_2000_py` (output)
9. `Automated_Grayscale.py` (output)
10. `2000plot.py` (output)

### Pipeline

The analysis pipeline involves several steps. To understand the flow of the pipeline I recommend:

1. To go through `source_file.py` (1), which loads constant values associated with the source and is called by other files.
2. Proceed to `effecient.py` (2), which is utilized by multiple other files.
3. `electrondensity.py` (3) calls `effecient.py`.
4. `emissivity.py` (4) calls `effecient.py`.
5. `freefreecalc.py` (5) is independent but is called by `abundance_calculator.py` (6).
6. `abundance_calculator.py` (6) serves as a semi-output file, interacting with most other files, and is called by output files.
7. `ionic_abundance_table.py` (7) produces output for specific sources, with the option for automation.
8. `GrayScale_2000_py` (8) produces semi-output for specific source transitions and types of grayscale. It is automated by `Automated_Grayscale.py` (9).
9. `2000plot.py` (10) is fully automated and generates plots for all available sources.

### Running the Analysis

To execute the analysis:

#### Abundance Plots
- For a specific transition of a specific source, use `GrayScale_2000_py` (8).
- For various sources, utilize `Automated_Grayscale.py` (9).
  read Note

#### Flux or Density Plots
- For specific transitions of specific sources, use `GrayScale_2000_py` (8).
- For various sources, utilize `Automated_Grayscale.py` (9).

#### Ionic Abundance Table
- For a specific source, modify the string value of `source` in line 20 of `ionic_abundance_table.py` (7) and adjust the saving path accordingly.
- For all sources, request modification.

#### Footprint Plots on Radio Images
- Use `2000plot.py` (10) and execute.

### Note
- For heavy calculations such as electron density and emissivity, pre-calculated values for W3 are available in `midspaxel.txt`.

---

Feel free to contact for questions!
