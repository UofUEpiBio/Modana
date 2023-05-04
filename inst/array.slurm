#!/bin/sh
#SBATCH --job-name=slurmRtry
#SBATCH --time=05:10:00 
# SBATCH --account=notchpeak-shared-short 
# SBATCH --partition=notchpeak-shared-short
# SBATCH --account=notchpeak-shared-short 
# SBATCH --partition=notchpeak-shared-short
# SBATCH --account=biostats-kp
# SBATCH --partition=biostats-kp
# SBATCH --account=greene
# SBATCH --partition=kingspeak
# SBATCH --account=greene
# SBATCH --partition=notchpeak
#SBATCH --account=greene
#SBATCH --partition=lonepeak
# SBATCH --account=smithp-guest
# SBATCH --partition=ash-guest
# SBATCH --out=finalarray
#SBATCH --err=finalarray
module load R/4.2.2
Rscript --vanilla modana.slurm.R
