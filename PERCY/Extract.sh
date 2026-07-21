#!/bin/bash
#SBATCH --mail-user=agogo-mawuli.percy@mayo.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=gene_extract_analysis
#SBATCH --partition=lg-n64-256g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --time=1-00:00:00
#SBATCH --mem=250G
#SBATCH --output=logs/%x.%N.%j.stdout
#SBATCH --error=logs/%x.%N.%j.stderr

set -euo pipefail
module load plink2/v2.00a.3LM.2022.0503

plink2 --pfile /fslustre/labs/ext_agogo_mawuli_percy_mayo_edu/SARP/SARP_FINAL/merged_qc_unrel\
       --extract range candidate_SNP_regions_1Mb.bed \
       --make-pgen \
       --out extracted_genes

