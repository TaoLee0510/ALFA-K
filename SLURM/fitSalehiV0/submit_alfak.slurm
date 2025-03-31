#! /bin/bash
mkdir -p logs

# Print how MYFILES is expanding before submitting the job
echo "Before submission - MYFILES contains:"
printf '%s\n' /home/4473331/projects/ALFA-K/data/salehi/alfak_inputs/*.Rds  # Debug raw glob expansion
echo "MYFILES expanded to: ${MYFILES[@]}"  # Debug: print list as expanded in bash

#efine MYFILES without modification (to see if the original is the problem)
MYFILES=(/home/4473331/projects/ALFA-K/data/salehi/alfak_inputs/*.Rds)
MYCOMMAND="Rscript /home/4473331/projects/ALFA-K/SLURM/fitSalehiV0/run_alfak.R"
MYMODULE="R"
NUMF="${#MYFILES[@]}"
SCRIPT=${0##*/} && SCRIPT=${SCRIPT%%.*}
CONC=${1:+%$1}

# Print debug info before submission
echo "SCRIPT: ${SCRIPT}"
echo "NUMF: ${NUMF} (number of files)"
echo "CONCURRENCY LIMIT: ${CONC}"
echo "COMMAND TO RUN: ${MYCOMMAND}"

#bmit to Slurm
sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1

#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-6
#SBATCH --array=0-$((${NUMF}-1))${CONC}
#SBATCH --job-name=${SCRIPT}
#SBATCH --output=logs/${SCRIPT}-%A_%a.out
#SBATCH --error=logs/${SCRIPT}-%A_%a.err

#  Inside Slurm job - debug outputs
echo "Inside Slurm job - SLURM_ARRAY_TASK_ID: \$SLURM_ARRAY_TASK_ID"
echo "SLURM_JOB_ID: \$SLURM_JOB_ID"
echo "Current working directory: \$(pwd)"
echo "Listing all files in input directory:"
ls -lh /home/4473331/projects/ALFA-K/data/salehi/alfak_inputs/  # Check if files exist

# Re-instantiate array inside job and print its values
MYFILES=(${MYFILES[@]})  # Expand inside Slurm job
echo "Inside Slurm - MYFILES expanded to: \${MYFILES[@]}"
echo "Number of files in MYFILES: \${#MYFILES[@]}"

# Select file based on SLURM_ARRAY_TASK_ID
THISFILE="\${MYFILES[\$SLURM_ARRAY_TASK_ID]}"
echo "Processing file: \${THISFILE}"

# Load R module and check if Rscript is available
ml ${MYMODULE}
which Rscript  # Verify Rscript is available

# Run the actual command but with debug wrapper
echo "Running: ${MYCOMMAND} \${THISFILE}"
${MYCOMMAND} "\${THISFILE}"

EOSUBMIT
