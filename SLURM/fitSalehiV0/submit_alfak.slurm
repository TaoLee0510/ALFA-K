#!/bin/bash

########################
# 1) Define your input files and concurrency
########################
MYFILES=(data/salehi/alfak_inputs/*.Rds)      # All .Rds files in alfak_inputs
MYCOMMAND="Rscript run_alfak.R"               # The command you want to run
MYMODULE="R/4.2.2"                            # Example module; adjust as needed
CONC=10                                       # Maximum number of parallel tasks
NUMF="${#MYFILES[@]}"                         # Count of input files

########################
# 2) Create a temporary SLURM submission (here-doc)
########################
sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --job-name=alfak
#SBATCH --output=logs/alfak_%A_%a.out
#SBATCH --error=logs/alfak_%A_%a.err
#SBATCH --array=0-$((${NUMF}-1))%${CONC}  # zero-based array, up to NUMF-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --partition=red       # or any other partition you prefer

########################
# 3) Load modules and run
########################
ml ${MYMODULE}               # Load R environment

MYFILES=(${MYFILES[@]})      # Re-instantiate array inside job
THISFILE="\${MYFILES[\$SLURM_ARRAY_TASK_ID]}"
echo "Running \${MYCOMMAND} on: \${THISFILE}"
echo "SLURM_ARRAY_TASK_ID: \$SLURM_ARRAY_TASK_ID"

# Actually run it
\${MYCOMMAND} "\${THISFILE}"

EOSUBMIT
