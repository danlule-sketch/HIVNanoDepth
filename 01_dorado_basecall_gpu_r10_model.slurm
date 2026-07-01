#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=nanoHIVDR01
#SBATCH --mail-user=daniel.lule@lshtm.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --output=dorado_rebasecall_01.out
#SBATCH --gres=gpu
#SBATCH --partition=gpu


# load the module, we use this for LSHTM HPC
module load dorado

cd /home/lsh2202849/

dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.2.0 FBD93406_daca3d03_253109ee_17.pod5 --kit-name SQK-NBD114-96 > /home/lsh2202849/FBD93406_daca3d03_253109ee_17.bam
