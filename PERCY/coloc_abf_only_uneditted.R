#!/usr/bin/env Rscript

#========================
# Colocalization pipeline: eQTL vs GWAS summary statistics
# Method: coloc.abf ONLY (Approximate Bayes Factor colocalization)
# No LD matrix required - uses summary statistics only
# Adapted for Google Cloud Platform (GCP)
# Input files:
#   - matrix_eqtl_output.txt : eQTL summary statistics
#   - ld_matched_gwas.txt    : GWAS summary statistics
#========================

#========================
# SECTION 1: GCP Package Installation
# Run once on a fresh GCP VM or Vertex AI notebook
#========================

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing("data.table")
install_if_missing("stringr")
install_if_missing("dplyr")
install_if_missing("yaml")
install_if_missing("BiocManager")
install_if_missing("snpStats", bioc = TRUE)
install_if_missing("coloc")

#========================
# SECTION 2: Load Libraries
#========================

suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(coloc))
suppressMessages(library(snpStats))
suppressMessages(library(dplyr))
suppressMessages(library(yaml))

options(scipen = 10)
options(datatable.fread.datatable = FALSE)

#========================
# SECTION 3: Global Parameters
#========================

THR_POST_PROB <- 0.75      # PP.H4 threshold for colocalization
GWAS_WINDOW   <- 500000    # 500 Kb window on each side of lead GWAS SNP
DUMMY_SE_VAL  <- 0.00001   # Replacement for zero or invalid SE entries
P12_PRIOR     <- 1e-5      # Prior probability of shared causal variant (coloc.abf)

# MHC region on chr6 - excluded from colocalization
MHC_chr6_LB <- 28477897
MHC_chr6_UB <- 33448354

#========================
# SECTION 4: File Paths
# Edit WORKING_DIR to match your GCP bucket mount or VM working directory
# e.g. "/mnt/gcs/your-bucket/Colocalization" for gcsfuse-mounted buckets
#========================

WORKING_DIR <- "."

eQTL_FILE <- file.path(WORKING_DIR, "matrix_eqtl_output.txt")
GWAS_FILE <- file.path(WORKING_DIR, "ld_matched_gwas.txt")
OUT_DIR   <- file.path(WORKING_DIR, "coloc_abf_only_output")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(OUT_DIR, "coloc_abf_run_summary.log")
sink(log_file, split = TRUE)

cat("=====================================================\n")
cat(" coloc.abf-Only Colocalization Pipeline - GCP Run\n")
cat(sprintf(" eQTL file  : %s\n", eQTL_FILE))
cat(sprintf(" GWAS file  : %s\n", GWAS_FILE))
cat(sprintf(" p12 prior  : %s\n", P12_PRIOR))
cat(sprintf(" Output dir : %s\n", OUT_DIR))
cat("=====================================================\n\n")

#========================
# SECTION 5: Read and Parse eQTL Data
# Expected columns in matrix_eqtl_output.txt:
#   ID, rsnumber, chr, chrom, pos, ref, alt, Ens_gene_id,
#   statistic, pvalue, FDR, beta, ALT_FREQS, MAF
#========================

cat("Reading eQTL data...\n")
eQTLData <- data.table::fread(eQTL_FILE, header = TRUE)
cat(sprintf("  Total eQTL entries loaded: %s\n", nrow(eQTLData)))

eQTLData <- eQTLData %>%
  rename(
    snp_id    = ID,
    rsID      = rsnumber,
    chr       = chr,
    chrom     = chrom,
    pos       = pos,
    ref       = ref,
    alt       = alt,
    eGeneID   = Ens_gene_id,
    statistic = statistic,
    pval_eQTL = pvalue,
    FDR_eQTL  = FDR,
    beta_eQTL = beta,
    AF_eQTL   = ALT_FREQS,
    MAF_eQTL  = MAF
  )

# Compute SE from beta and two-tailed p-value z-score approximation
# Reference: https://www.biostars.org/p/431875/
eQTLData$SE_eQTL <- abs(eQTLData$beta_eQTL / qnorm(eQTLData$pval_eQTL / 2))

idx_bad_se <- which(eQTLData$SE_eQTL == 0 |
                    is.nan(eQTLData$SE_eQTL) |
                    is.infinite(eQTLData$SE_eQTL))
if (length(idx_bad_se) > 0) {
  eQTLData$SE_eQTL[idx_bad_se] <- DUMMY_SE_VAL
}

if (!grepl("chr", eQTLData$chr[1])) {
  eQTLData$chr <- paste0("chr", eQTLData$chr)
}

cat(sprintf("  eQTL entries after preprocessing: %s\n", nrow(eQTLData)))

#========================
# SECTION 6: Read and Parse GWAS Data
# Expected format of ld_matched_gwas.txt (no header, space-delimited):
#   Col1=SNP_ID, Col2=chr_num, Col3=pos, Col4=alt, Col5=ref,
#   Col6=beta, Col7=SE, Col8=z_stat, Col9=pvalue
#========================

cat("\nReading GWAS data...\n")
gwas_raw <- data.table::fread(GWAS_FILE, header = FALSE)
cat(sprintf("  GWAS columns detected: %s\n", ncol(gwas_raw)))
cat(sprintf("  GWAS rows loaded     : %s\n", nrow(gwas_raw)))

colnames(gwas_raw) <- c("snp_id", "chrom_num", "pos", "alt", "ref",
                         "beta_GWAS", "SE_gwas", "z_stat", "pval_GWAS")

GWASData <- gwas_raw %>%
  mutate(
    chr       = paste0("chr", chrom_num),
    pos       = as.integer(pos),
    beta_GWAS = as.numeric(beta_GWAS),
    SE_gwas   = as.numeric(SE_gwas),
    pval_GWAS = as.numeric(pval_GWAS)
  )

idx_gwas_bad <- which(GWASData$SE_gwas == 0 | is.na(GWASData$SE_gwas))
if (length(idx_gwas_bad) > 0) {
  GWASData$SE_gwas[idx_gwas_bad] <- DUMMY_SE_VAL
}

cat(sprintf("  GWAS entries after preprocessing: %s\n", nrow(GWASData)))

#========================
# SECTION 7: Define GWAS Loci
# Clump GWAS SNPs at p < 5e-8; define 500 Kb windows around each lead SNP
#========================

ColocSNPInfoFile          <- file.path(OUT_DIR, "FINAL_coloc_abf_Summary.txt")
ColocSNPInfoFile_credible <- file.path(OUT_DIR, "FINAL_coloc_abf_95pct_credible_set.txt")
ColocFailedLog            <- file.path(OUT_DIR, "coloc_abf_Failed_Loci_Genes.txt")

bool_Coloc_Summary_DF    <- FALSE
bool_Credible_Summary_DF <- FALSE
failed_entries           <- list()

GWASChrList <- as.vector(sort(unique(GWASData$chr)))

for (chridx in seq_along(GWASChrList)) {

  currchr <- GWASChrList[chridx]
  cat(sprintf("\n\n========== Processing chromosome: %s ==========\n", currchr))

  GWASData_currchr <- GWASData[GWASData$chr == currchr, ]
  eQTLData_currchr <- eQTLData[eQTLData$chr == currchr, ]
  cat(sprintf("  GWAS SNPs : %s\n", nrow(GWASData_currchr)))
  cat(sprintf("  eQTL entries : %s\n", nrow(eQTLData_currchr)))

  if (nrow(GWASData_currchr) == 0 || nrow(eQTLData_currchr) == 0) {
    cat("  No data on this chromosome. Skipping.\n")
    next
  }

  # Clump GWAS loci
  gwasdata_sorted <- GWASData_currchr[order(GWASData_currchr$pval_GWAS), ]
  bool_GWAS_Loci  <- FALSE

  while (nrow(gwasdata_sorted) > 0) {
    if (gwasdata_sorted$pval_GWAS[1] > 5e-8) break

    lead_pos <- gwasdata_sorted$pos[1]
    cat(sprintf("  Lead GWAS SNP: %s  pos: %s\n", currchr, lead_pos))

    startpos   <- max(0, lead_pos - GWAS_WINDOW)
    endpos     <- lead_pos + GWAS_WINDOW
    currLociDF <- data.frame(chr = currchr, start = startpos, end = endpos)

    if (!bool_GWAS_Loci) {
      GWAS_Loci_DF <- currLociDF
      bool_GWAS_Loci <- TRUE
    } else {
      GWAS_Loci_DF <- rbind(GWAS_Loci_DF, currLociDF)
    }

    idx_remove <- which(gwasdata_sorted$pos >= startpos &
                        gwasdata_sorted$pos <= endpos)
    gwasdata_sorted <- gwasdata_sorted[-idx_remove, ]
  }

  if (!bool_GWAS_Loci) {
    cat("  No genome-wide significant loci on this chromosome. Skipping.\n")
    next
  }

  cat(sprintf("  GWAS loci defined: %s\n", nrow(GWAS_Loci_DF)))

  #========================
  # SECTION 8: Process Each GWAS Locus
  #========================

  for (lociidx in seq_len(nrow(GWAS_Loci_DF))) {

    startpos    <- GWAS_Loci_DF$start[lociidx]
    endpos      <- GWAS_Loci_DF$end[lociidx]
    locus_label <- paste0(currchr, "_", startpos, "_", endpos)
    cat(sprintf("\n  --- Locus: %s ---\n", locus_label))

    currloci_GWASdata <- GWASData_currchr[
      GWASData_currchr$pos >= startpos & GWASData_currchr$pos <= endpos, ]
    currloci_eqtldata <- eQTLData_currchr[
      eQTLData_currchr$pos >= startpos & eQTLData_currchr$pos <= endpos, ]

    cat(sprintf("    GWAS SNPs in locus : %s\n", nrow(currloci_GWASdata)))
    cat(sprintf("    eQTL entries in locus : %s\n", nrow(currloci_eqtldata)))

    if (nrow(currloci_GWASdata) == 0 || nrow(currloci_eqtldata) == 0) {
      cat("    Empty locus after subsetting. Skipping.\n")
      next
    }

    # Merge eQTL and GWAS on chr + pos
    merge_data <- dplyr::inner_join(currloci_eqtldata, currloci_GWASdata,
                                    by = c("chr", "pos"))

    if (nrow(merge_data) == 0) {
      cat("    No overlapping SNPs between eQTL and GWAS. Skipping.\n")
      next
    }

    # Exclude MHC region
    if (currchr == "chr6") {
      mhc_idx <- which(merge_data$pos >= MHC_chr6_LB &
                       merge_data$pos <= MHC_chr6_UB)
      if (length(mhc_idx) > 0) merge_data <- merge_data[-mhc_idx, ]
    }

    # Remove rows with NA in key fields
    na_idx <- which(
      is.na(merge_data$pval_GWAS)  | is.na(merge_data$pval_eQTL) |
      is.na(merge_data$beta_GWAS)  | is.na(merge_data$SE_gwas)   |
      is.na(merge_data$beta_eQTL) | is.na(merge_data$SE_eQTL)
    )
    if (length(na_idx) > 0) merge_data <- merge_data[-na_idx, ]

    if (nrow(merge_data) == 0) {
      cat("    No valid entries after NA filtering. Skipping.\n")
      next
    }

    cat(sprintf("    Merged SNPs for colocalization: %s\n", nrow(merge_data)))

    # Create per-locus output directory
    CurrLoci_OutDir <- file.path(OUT_DIR, paste0("Locus_", locus_label))
    dir.create(CurrLoci_OutDir, showWarnings = FALSE, recursive = TRUE)

    write.table(merge_data,
                file.path(CurrLoci_OutDir, "merged_GWAS_eQTL_input.txt"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    #========================
    # SECTION 9: Process Each Gene in This Locus
    #========================

    gene_list <- unique(merge_data$eGeneID)
    cat(sprintf("    Genes to process: %s\n", length(gene_list)))

    for (gene in gene_list) {

      gene_data <- merge_data[merge_data$eGeneID == gene, ]
      cat(sprintf("\n    Gene: %s  SNPs: %s\n", gene, nrow(gene_data)))

      if (nrow(gene_data) < 2) {
        cat("    Too few SNPs. Skipping.\n")
        failed_entries <- c(failed_entries,
                            list(data.frame(locus  = locus_label,
                                            gene   = gene,
                                            reason = "Too few SNPs (< 2)")))
        next
      }

      n_gwas <- nrow(currloci_GWASdata)
      n_eqtl <- nrow(currloci_eqtldata[currloci_eqtldata$eGeneID == gene, ])

      #========================
      # SECTION 10: Build coloc Datasets
      #========================

      dataset_gwas <- list(
        beta    = gene_data$beta_GWAS,
        varbeta = gene_data$SE_gwas ^ 2,
        pvalues = gene_data$pval_GWAS,
        type    = "quant",
        N       = max(n_gwas, 1),
        sdY     = 1
      )

      dataset_eqtl <- list(
        beta    = gene_data$beta_eQTL,
        varbeta = gene_data$SE_eQTL ^ 2,
        pvalues = gene_data$pval_eQTL,
        type    = "quant",
        N       = max(n_eqtl, 1),
        MAF     = gene_data$MAF_eQTL,
        sdY     = 1
      )

      #========================
      # SECTION 11: Run coloc.abf
      #========================

      gene_out_dir <- file.path(CurrLoci_OutDir, gene)
      dir.create(gene_out_dir, showWarnings = FALSE, recursive = TRUE)

      coloc_result <- NULL
      cat("    Running coloc.abf...\n")

      tryCatch({
        coloc_result <- coloc::coloc.abf(
          dataset1 = dataset_gwas,
          dataset2 = dataset_eqtl,
          MAF      = gene_data$MAF_eQTL,
          p12      = P12_PRIOR
        )
      }, error = function(e) {
        cat(sprintf("    coloc.abf failed: %s\n", e$message))
        failed_entries <<- c(failed_entries,
                             list(data.frame(locus  = locus_label,
                                             gene   = gene,
                                             reason = paste("coloc.abf failed:", e$message))))
      })

      if (is.null(coloc_result)) {
        cat("    Skipping gene (coloc.abf failed).\n")
        next
      }

      cat("    coloc.abf completed successfully.\n")

      #========================
      # SECTION 12: Save Per-Gene coloc.abf Results
      #========================

      # Save summary posterior probabilities
      write.table(coloc_result$summary,
                  file.path(gene_out_dir, "coloc_abf_summary.txt"),
                  row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

      # Save per-SNP posterior probabilities
      write.table(coloc_result$results,
                  file.path(gene_out_dir, "coloc_abf_results_per_snp.txt"),
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

      pp_H0 <- as.numeric(coloc_result$summary["PP.H0.abf"])
      pp_H1 <- as.numeric(coloc_result$summary["PP.H1.abf"])
      pp_H2 <- as.numeric(coloc_result$summary["PP.H2.abf"])
      pp_H3 <- as.numeric(coloc_result$summary["PP.H3.abf"])
      pp_H4 <- as.numeric(coloc_result$summary["PP.H4.abf"])

      cat(sprintf("    PP.H0=%.4f  PP.H1=%.4f  PP.H2=%.4f  PP.H3=%.4f  PP.H4=%.4f\n",
                  ifelse(is.na(pp_H0),0,pp_H0),
                  ifelse(is.na(pp_H1),0,pp_H1),
                  ifelse(is.na(pp_H2),0,pp_H2),
                  ifelse(is.na(pp_H3),0,pp_H3),
                  ifelse(is.na(pp_H4),0,pp_H4)))

      #========================
      # SECTION 13: Extract 95% Credible Set and Collect Signals
      #========================

      if (!is.na(pp_H4) && pp_H4 > THR_POST_PROB) {
        cat(sprintf("    COLOCALIZATION FOUND for gene %s (PP.H4 = %.4f)\n", gene, pp_H4))

        res_df    <- coloc_result$results
        res_df    <- res_df[order(res_df$SNP.PP.H4, decreasing = TRUE), ]
        cs_cumsum <- cumsum(res_df$SNP.PP.H4)
        w         <- which(cs_cumsum > 0.95)[1]
        if (is.na(w)) w <- nrow(res_df)
        credible_set <- res_df[1:w, ]

        # Tag with gene and locus metadata
        credible_set$gene  <- gene
        credible_set$locus <- locus_label
        credible_set$PP_H4 <- pp_H4

        write.table(credible_set,
                    file.path(gene_out_dir, "coloc_abf_95pct_credible_set.txt"),
                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

        # Top colocalized SNP summary row
        top_snp_row          <- credible_set[1, , drop = FALSE]
        top_snp_row$gene     <- gene
        top_snp_row$locus    <- locus_label
        top_snp_row$PP_H0    <- pp_H0
        top_snp_row$PP_H1    <- pp_H1
        top_snp_row$PP_H2    <- pp_H2
        top_snp_row$PP_H3    <- pp_H3
        top_snp_row$PP_H4    <- pp_H4
        top_snp_row$n_credible_snps <- nrow(credible_set)

        if (!bool_Coloc_Summary_DF) {
          Final_topColocDF      <- top_snp_row
          bool_Coloc_Summary_DF <- TRUE
        } else {
          Final_topColocDF <- rbind(Final_topColocDF, top_snp_row)
        }

        if (!bool_Credible_Summary_DF) {
          Final_credibleColocDF      <- credible_set
          bool_Credible_Summary_DF   <- TRUE
        } else {
          Final_credibleColocDF <- rbind(Final_credibleColocDF, credible_set)
        }

      } else {
        cat(sprintf("    No colocalization signal for gene %s (PP.H4 = %.4f)\n",
                    gene, ifelse(is.na(pp_H4), 0, pp_H4)))
      }

    }  # end gene loop

  }  # end locus loop

}  # end chromosome loop

#========================
# SECTION 14: Write Final Summary Files
#========================

cat("\n\n========== Writing final output files ==========\n")

if (bool_Coloc_Summary_DF && exists("Final_topColocDF") && nrow(Final_topColocDF) > 0) {
  write.table(Final_topColocDF, ColocSNPInfoFile,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  cat(sprintf("  Top colocalized SNP summary written: %s\n", ColocSNPInfoFile))
} else {
  cat("  No colocalization signals found above threshold.\n")
}

if (bool_Credible_Summary_DF && exists("Final_credibleColocDF") && nrow(Final_credibleColocDF) > 0) {
  write.table(Final_credibleColocDF, ColocSNPInfoFile_credible,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  cat(sprintf("  95pct credible set summary written: %s\n", ColocSNPInfoFile_credible))
}

if (length(failed_entries) > 0) {
  failed_df <- do.call(rbind, failed_entries)
  write.table(failed_df, ColocFailedLog,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  cat(sprintf("  Failed loci/genes log written: %s\n", ColocFailedLog))
}

cat("\n coloc.abf-only pipeline complete.\n")
sink()
