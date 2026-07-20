#!/bin/bash
#SBATCH --mail-user=agogo-mawuli.percy@mayo.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=LD_matrix
#SBATCH --partition=med-n16-64g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --mem=50G
#SBATCH --output=logs/%x.%N.%j.stdout
#SBATCH --error=logs/%x.%N.%j.stderr



set -euo pipefail

source /research/bsi/tools/biotools/conda/24.3.0/etc/profile.d/conda.sh
conda activate genius

which plink2
plink2 --version

mkdir -p /fslustre/labs/ext_agogo_mawuli_percy_mayo_edu/SARP/SARP_FINAL/


plink2 \
 --pfile /fslustre/labs/ext_agogo_mawuli_percy_mayo_edu/SARP/extracted_genes \
 --keep-allele-order \
 --r-unphased square \
 --out locus_SARP
