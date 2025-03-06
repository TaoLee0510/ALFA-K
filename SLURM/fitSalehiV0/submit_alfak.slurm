#! /bin/bash
mkdir -p logs
MYFILES=(/home/4473331/projects/ALFA-K/data/salehi/alfak_inputs/*.Rds) # array of files in dir
MYCOMMAND="Rscript /home/4473331/projects/ALFA-K/SLURM/fitSalehiV0/run_alfak.R"                    # run this command
MYMODULE="R"                                  # use this software module

# you might not have to change these
NUMF="${#MYFILES[@]}"                             # filecount in array
SCRIPT=${0##*/} && SCRIPT=${SCRIPT%%.*}           # remove path and ext
CONC=${1:+%$1}                                    # Concurrency cmd line

sbatch <<EOSUBMIT
#! /bin/bash
#SBATCH --nodes 1                       # one node
#SBATCH --ntasks 1                      # one task
#SBATCH --cpus-per-task 50               # one cpu core per task
#SBATCH --mem-per-cpu 4G                # 4G of memory per core
#SBATCH --time 1-6                      # kill job after 1 day and 6 hours
#SBATCH --array 0-$((${NUMF}-1))${CONC} # zero based array jobs with concurrency
#SBATCH --job-name ${SCRIPT}            # Jobname based on script name
#SBATCH --output=logs/${SCRIPT}-%A_%a.out    # Store STDOUT in logs/ directory
#SBATCH --error=logs/${SCRIPT}-%A_%a.err     # Store STDERR in logs/ directory


MYFILES=(${MYFILES[@]})                 # re-instantiate array of data files
ml ${MYMODULE}                          # load our environment module
which ${MYCOMMAND} 2> /dev/null         # print path/version of our software

echo "${MYCOMMAND} \${MYFILES[\$SLURM_ARRAY_TASK_ID]}"

echo "NUMF: ${NUMF}"
echo "SCRIPT: ${SCRIPT}"
echo "SLURM_JOB_ID: \${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: \${SLURM_ARRAY_TASK_ID}"
echo "SLURM_ARRAY_TASK_MAX: \${SLURM_ARRAY_TASK_MAX}"
echo "SLURM_ARRAY_TASK_MIN: \${SLURM_ARRAY_TASK_MIN}"
echo "SLURM_ARRAY_TASK_COUNT: \${SLURM_ARRAY_TASK_COUNT}"
echo "SLURM_SUBMIT_DIR: \${SLURM_SUBMIT_DIR}"
echo "SLURMD_NODENAME: \${SLURMD_NODENAME}"
echo "SLURM_JOB_CPUS_PER_NODE: \${SLURM_JOB_CPUS_PER_NODE}"
echo "SLURM_JOB_NODELIST: \${SLURM_JOB_NODELIST}"
echo "SLURM_MEM_PER_CPU: \${SLURM_MEM_PER_CPU}"

EOSUBMIT