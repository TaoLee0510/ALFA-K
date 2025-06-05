# Chromosome Fitness Landscapes – alfak Analysis Code
ALFA-K takes longitudinal single cell sequencing data from an evolving
cell population as input, then estimates a local fitness landscape
encompassing thousands of karyotypes located near to the input data in
karyotype space. This repository contains source code and examples for
running ALFA-K, as well as an Agent Based Model (ABM) which simulates
evolving cell populations using fitness landscapes estimated by ALFA-K.
This repository contains all scripts and data needed to reproduce the analyses and figures from:

Beck, Richard J., and Noemi Andor. “Local Adaptive Mapping of Karyotype
Fitness Landscapes.” bioRxiv (2023): 2023-07. (pending citation update)

---

Quick Start
-----------

Clone the repository and recreate the exact R environment:

    git clone https://github.com/Richard-Beck/ALFA-K.git
    cd ALFA-K

    # Set up R environment
    install.packages("renv")
    renv::restore()

    # Run the pipeline
    Rscript scripts/01_preprocess_cell_lines.R
    # See note about batch scripts below
    ./scripts/02_batch_run_alfak.R
    # ... etc.

All outputs will be written to `data/` and `figs/`.

---

Repository Structure
--------------------

.
├── README.md             # You're here  
├── install.R             # Calls renv::restore()  
├── renv.lock             # Exact package versions  
├── scripts/              # Main analysis scripts  
├── R/                    # Utility functions (themes, packages, etc.)  
├── data/  
│   ├── raw/              # Immutable input data  
│   └── processed/        # Outputs from scripts   
├── figures/              # Generated figures 


---

Analysis Details
-----------------

Scripts should be run in order. 

Some scripts are standalone, others 


---

Software Environment
--------------------

- R ≥ 4.3  
- alfakR (https://github.com/Richard-Beck/alfakR, ≥ 0.9.1)  
- All required packages and versions are listed in `renv.lock`





  



