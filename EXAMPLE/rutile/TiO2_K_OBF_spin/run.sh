#!/bin/bash
#SBATCH --partition=nano
#SBATCH --account=nano
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --job-name=tio2qe
#SBATCH -e stderr
#SBATCH -o stdout
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pcbee912@gmail.com

cd $SLURM_SUBMIT_DIR

/global/home/users/yufengl/ocean_bin/ocean.pl ./rutile.in > rutile.out


