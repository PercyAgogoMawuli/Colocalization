#!/bin/bash

awk ' FNR==NR {
    # 1. Load unique LD matrix SNPs
    ld_snps[$1]; next
}
NR==1 {
    # 2. Print header
    print $0; next
}
($1 in ld_snps) {
    # 3. Only print if we havent seen this SNP in the GWAS file yet
    if (!seen[$1]++) { print $0
    }
}' locus_SARP.unphased.vcor1.vars formatted_gwas.txt > ld_matched_gwas.txt
