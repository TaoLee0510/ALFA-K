#!/bin/bash
#SBATCH --job-name=alfak
#SBATCH --output=logs/alfak_%A_%a.out
#SBATCH --error=logs/alfak_%A_%a.err
#SBATCH --array=1-1  # Placeholder, overridden at submission
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50  # Match num_cores inside alfak()
#SBATCH --mem=200G          # Adjust based on actual memory usage
#SBATCH --time=12:00:00    # Set appropriate runtime

# Load R environment (adjust module as needed)
module load R

# Run R script for this task
Rscript run_alfak.R
