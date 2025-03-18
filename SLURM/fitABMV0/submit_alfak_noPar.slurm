#! /bin/bash
mkdir -p logs

# Define constants
MYCOMMAND="Rscript /home/4473331/projects/ALFA-K/SLURM/fitABMV0/run_alfak_noPar.R"
MYMODULE="R"
NUMF=400  # Set to 400 instead of counting files
SCRIPT=${0##*/} && SCRIPT=${SCRIPT%%.*}
CONC=${1:+%$1}

# Print debug info before submission
echo "SCRIPT: ${SCRIPT}"
echo "NUMF: ${NUMF} (number of tasks)"
echo "CONCURRENCY LIMIT: ${CONC}"
echo "COMMAND TO RUN: ${MYCOMMAND}"

# Submit to Slurm
sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-6
#SBATCH --qos=small
#SBATCH --array=1-${NUMF}${CONC}
#SBATCH --job-name=${SCRIPT}
#SBATCH --output=logs/${SCRIPT}-%A_%a.out
#SBATCH --error=logs/${SCRIPT}-%A_%a.err

# Inside Slurm job - debug outputs
echo "Inside Slurm job - SLURM_ARRAY_TASK_ID: \$SLURM_ARRAY_TASK_ID"
echo "SLURM_JOB_ID: \$SLURM_JOB_ID"
echo "Current working directory: \$(pwd)"

# Load R module and check if Rscript is available
ml ${MYMODULE}
which Rscript  # Verify Rscript is available

# Run the actual command but with debug wrapper
echo "Running: ${MYCOMMAND} \$SLURM_ARRAY_TASK_ID"
${MYCOMMAND} "\$SLURM_ARRAY_TASK_ID"

EOSUBMIT
