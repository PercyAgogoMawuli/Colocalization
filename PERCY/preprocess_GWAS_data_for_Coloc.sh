#!/bin/bash


awk ' BEGIN { printf "SNP\tCHR\tPOS\tA1\tA2\tOR\tSE\tSTAT\tP\n" } NR==1 { for(i=1; i<=NF; i++) { if($i=="#CHROM") chr_idx=i; if($i=="POS") pos_idx=i; if($i=="ID") snp_idx=i; if($i=="A1") a1_idx=i; if($i=="REF") a2_idx=i; 
        if($i=="OR") or_idx=i; if($i=="LOG(OR)_SE") se_idx=i; if($i=="Z_STAT") stat_idx=i; if($i=="P") p_idx=i;
    }
    next;
}
{ print $snp_idx, $chr_idx, $pos_idx, $a1_idx, $a2_idx, $or_idx, $se_idx, $stat_idx, $p_idx;
}' OFS="\t" gwas_clean.txt > formatted_gwas.txt
