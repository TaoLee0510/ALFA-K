#! /bin/bash
#SBATCH --ntasks 1            # 1 task, gets 1 node and 1 cpu by default
#SBATCH --time 0-2:30:00      # kill job after 0 days, 2h and 30 mins
#SBATCH --output myjob-%j.out # write to custom output myjob-xxx.out
ml R                      # select the version of R to be used
which R                       # print the R installation/version used here 
Rscript run_alfak.R /home/4473331/projects/ALFA-K/data/salehi/alfak_inputs/SA906a_X50_l_5_d1_1_d2_0.Rds             # launch the script 
