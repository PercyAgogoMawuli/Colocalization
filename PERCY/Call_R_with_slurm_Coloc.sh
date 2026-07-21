#!/bin/bash
#SBATCH --mail-user=agogo-mawuli.percy@mayo.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=R
#SBATCH --partition=huge-n128-512g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=115
#SBATCH --time=7-00:00:00
#SBATCH --mem=400G
#SBATCH --output=logs/%x.%N.%j.stdout
#SBATCH --error=logs/%x.%N.%j.stderr

set -euo pipefail
module load r/4.5.1

Rscript /fslustre/labs/ext_agogo_mawuli_percy_mayo_edu/SARP/SARP_FINAL/Colocalization/candidate_coloc_susie_only.R


