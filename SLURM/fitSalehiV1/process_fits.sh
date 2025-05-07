#! /bin/bash
#submit to Slurm
sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1

#SBATCH --cpus-per-task=51
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-6
#SBATCH --qos=small


# Load R module and check if Rscript is available
ml R
which Rscript  # Verify Rscript is available

Rscript process_salehi_sweep_data.R

EOSUBMIT


