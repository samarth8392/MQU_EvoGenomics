#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=220G
#SBATCH -t 14-00:00:00
#SBATCH --job-name=1k.expand2
#SBATCH -e 1k.expand2
#SBATCH -o /scratch/bell/mathur20/MQU/ch3/revise/slim/1k/output/1k.expand2

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/04/20                 Last Modified: 02/15/22 	###
###########################################################################
###########################################################################
###                     submit_slim_jobs.sh                   			###
###########################################################################

cd $SLURM_SUBMIT_DIR


cd /scratch/bell/mathur20/MQU/ch3/revise/slim/1k
~/build/slim sims/1k_revise.R3.expand.slim
