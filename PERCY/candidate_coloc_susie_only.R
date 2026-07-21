#!/usr/bin/env Rscript

#========================
# Colocalization pipeline: eQTL vs GWAS summary statistics
# Method: SuSiE ONLY (coloc.susie with LD matrix)
# No fallback to coloc.abf - if SuSiE fails, the locus/gene is skipped and logged
# Adapted for Google Cloud Platform (GCP)
# Input files:
#   - matrix_eqtl_output.txt         : eQTL summary statistics
#   - ld_matched_gwas.txt             : GWAS summary statistics (LD-matched)
#   - locus_SARP.unphased.vcor1       : LD correlation matrix
#   - locus_SARP.unphased.vcor1.vars  : SNP IDs for LD matrix rows/columns
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
install_if_missing("Matrix")
install_if_missing("susieR")
install_if_missing("BiocManager")
install_if_missing("snpStats", bioc = TRUE)
install_if_missing("coloc")

#========================
# SECTION 2: Load Libraries
#========================
# ----- Set working directory -----
setwd("/fslustre/labs/ext_agogo_mawuli_percy_mayo_edu/SARP/SARP_FINAL/Colocalization")  # <-- adjust as needed
cat("Working directory set to:", getwd(), "\n")


suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(coloc))
suppressMessages(library(snpStats))
suppressMessages(library(dplyr))
suppressMessages(library(yaml))
suppressMessages(library(Matrix))
suppressMessages(library(susieR))

options(scipen = 10)
options(datatable.fread.datatable = FALSE)

#========================
# SECTION 3: Global Parameters
#========================

THR_POST_PROB <- 0.75      # PP.H4 threshold for colocalization
GWAS_WINDOW   <- 1000000    # 1Mb  window on each side of lead GWAS SNP
DUMMY_SE_VAL  <- 0.00001   # Replacement for zero or invalid SE entries

# MHC region on chr6 - excluded from colocalization
MHC_chr6_LB <- 28477897
MHC_chr6_UB <- 33448354

#========================
# SECTION 4: File Paths
# Edit WORKING_DIR to match your GCP bucket mount or VM working directory
# e.g. "/mnt/gcs/your-bucket/Colocalization" for gcsfuse-mounted buckets
#========================

WORKING_DIR <- "/fslustre/labs/ext_agogo_mawuli_percy_mayo_edu/SARP/SARP_FINAL/Colocalization"

eQTL_FILE   <- file.path(WORKING_DIR, "matrix_eqtl_output.txt")
GWAS_FILE   <- file.path(WORKING_DIR, "ld_matched_gwas.txt")
LD_MAT_FILE <- file.path(WORKING_DIR, "locus_SARP.unphased.vcor1")
LD_VAR_FILE <- file.path(WORKING_DIR, "locus_SARP.unphased.vcor1.vars")
OUT_DIR     <- file.path(WORKING_DIR, "susie_only_output")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(OUT_DIR, "susie_run_summary.log")
sink(log_file, split = TRUE)

cat("=====================================================\n")
cat(" SuSiE-Only Colocalization Pipeline - GCP Run\n")
cat(sprintf(" eQTL file  : %s\n", eQTL_FILE))
cat(sprintf(" GWAS file  : %s\n", GWAS_FILE))
cat(sprintf(" LD matrix  : %s\n", LD_MAT_FILE))
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


#addition
#========================
# Read replicated SNP list
#========================

REP_SNP_FILE <- file.path(WORKING_DIR,
                          "replicated_snps.txt")

RepSNPs <- fread(REP_SNP_FILE,
                 header = TRUE)

if(!"chr" %in% names(RepSNPs)){
    RepSNPs$chr <- paste0("chr",
                          gsub("chr","",RepSNPs$chr))
}

#========================
# SECTION 7: Read LD Matrix
# locus_SARP.unphased.vcor1      : square correlation matrix, no header or rownames
# locus_SARP.unphased.vcor1.vars : one SNP ID per line (row/col labels for matrix)
#========================

cat("\nReading LD matrix...\n")
ld_vars <- trimws(readLines(LD_VAR_FILE))
cat(sprintf("  LD matrix SNPs: %s\n", length(ld_vars)))

ld_mat_raw <- data.table::fread(LD_MAT_FILE, header = FALSE)
LD_matrix  <- as.matrix(ld_mat_raw)
rownames(LD_matrix) <- ld_vars
colnames(LD_matrix) <- ld_vars
cat(sprintf("  LD matrix dimensions: %s x %s\n", nrow(LD_matrix), ncol(LD_matrix)))

# Symmetrize and enforce unit diagonal (numerical precision fix)
LD_matrix <- (LD_matrix + t(LD_matrix)) / 2
diag(LD_matrix) <- 1

#========================
# SECTION 8: Define loci from replicated SNPs
#========================

RepSNPs_chr <- RepSNPs[
    RepSNPs$chr == currchr,
]

if(nrow(RepSNPs_chr)==0){

    cat("  No replicated SNPs on this chromosome.\n")
    next
}

GWAS_Loci_DF <- data.frame()

for(i in seq_len(nrow(RepSNPs_chr))){

    lead_pos <- RepSNPs_chr$pos[i]

    startpos <- max(0,
                    lead_pos-GWAS_WINDOW)

    endpos <- lead_pos+GWAS_WINDOW

    GWAS_Loci_DF <- rbind(
        GWAS_Loci_DF,
        data.frame(

            chr=currchr,

            lead_snp=RepSNPs_chr$snp_id[i],

            lead_pos=lead_pos,

            start=startpos,

            end=endpos,

            stringsAsFactors=FALSE
        )
    )

}


cat(sprintf("  GWAS loci defined: %s\n", nrow(GWAS_Loci_DF)))




merge_loci <- function(loci){

    loci <- loci[order(loci$start), ]

    merged <- list()

    current <- loci[1, ]

    current$anchor_snps <- current$lead_snp

    for(i in 2:nrow(loci)){

        next_locus <- loci[i, ]

        if(next_locus$start <= current$end){

            current$end <- max(current$end,
                               next_locus$end)

            current$anchor_snps <- paste(
                current$anchor_snps,
                next_locus$lead_snp,
                sep=";"
            )

        } else{

            merged[[length(merged)+1]] <- current

            current <- next_locus

            current$anchor_snps <- current$lead_snp
        }
    }

    merged[[length(merged)+1]] <- current

    do.call(rbind, merged)

}

GWAS_Loci_DF <- merge_loci(GWAS_Loci_DF)



  #========================
  # SECTION 9: Process Each GWAS Locus
  #========================

  for (lociidx in seq_len(nrow(GWAS_Loci_DF))) {

    startpos    <- GWAS_Loci_DF$start[lociidx]
    endpos      <- GWAS_Loci_DF$end[lociidx]
    lead_snp <- GWAS_Loci_DF$lead_snp[lociidx]
    anchor_snps <- GWAS_Loci_DF$anchor_snps[lociidx]

    #locus_label <- paste0(
     #   lead_snp,
      #  "_",
       # startpos,
        #"_",
        #endpos
     #)
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

     ####################################################
     # Harmonize alleles between GWAS and eQTL
     ####################################################

     # Same allele orientation
     same_orientation <-
         merge_data$ref.x == merge_data$ref.y &
         merge_data$alt.x == merge_data$alt.y

     # Reverse allele orientation
     reverse_orientation <-
         merge_data$ref.x == merge_data$alt.y &
         merge_data$alt.x == merge_data$ref.y

     # Keep only SNPs with matching alleles
     merge_data <- merge_data[
         same_orientation | reverse_orientation,
     ]

     # Flip eQTL effect sizes when alleles are reversed
     merge_data$beta_eQTL[reverse_orientation] <-
         -merge_data$beta_eQTL[reverse_orientation]

     cat(sprintf(
         "    SNPs after allele harmonization: %s\n",
         nrow(merge_data)
     ))

     if (nrow(merge_data) == 0) {
         cat("    No allele-matched SNPs. Skipping.\n")
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
    # SECTION 10: Match LD Matrix to Merged SNPs
    #========================

    snp_col <- if ("snp_id.x" %in% colnames(merge_data)) "snp_id.x" else
               if ("snp_id"   %in% colnames(merge_data)) "snp_id"   else NULL

    if (is.null(snp_col)) {
      cat("    Cannot identify SNP ID column in merged data. Skipping locus.\n")
      failed_entries <- c(failed_entries,
                          list(data.frame(locus = locus_label, gene = "ALL",
                                          reason = "No SNP ID column found")))
      next
    }

    snp_ids_in_data <- merge_data[[snp_col]]
    shared_snps     <- intersect(snp_ids_in_data, ld_vars)
    cat(sprintf("    SNPs matching LD matrix: %s / %s\n",
                length(shared_snps), nrow(merge_data)))

    if (length(shared_snps) < 2) {
      cat("    Insufficient LD overlap (< 2 SNPs). Skipping locus.\n")
      failed_entries <- c(failed_entries,
                          list(data.frame(locus = locus_label, gene = "ALL",
                                          reason = "Insufficient LD overlap")))
      next
    }

    # Subset merge_data and LD matrix to shared SNPs
    merge_data <- merge_data[snp_ids_in_data %in% shared_snps, ]
    snp_order  <- merge_data[[snp_col]]
    LD_sub     <- LD_matrix[snp_order, snp_order, drop = FALSE]

    # Ensure positive semi-definiteness (required by SuSiE)
    eig_vals <- eigen(LD_sub, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig_vals < -1e-6)) {
      cat("    LD matrix not PSD - applying eigenvalue regularization.\n")
      eig_decomp <- eigen(LD_sub, symmetric = TRUE)
      eig_decomp$values[eig_decomp$values < 0] <- 0
      LD_sub <- eig_decomp$vectors %*%
                diag(eig_decomp$values) %*%
                t(eig_decomp$vectors)
      LD_sub <- (LD_sub + t(LD_sub)) / 2
      diag(LD_sub) <- 1
    }

    #========================
    # SECTION 11: Process Each Gene in This Locus
    #========================

    gene_list <- unique(merge_data$eGeneID)
    cat(sprintf("    Genes to process: %s\n", length(gene_list)))

    for (gene in gene_list) {

      gene_data <- merge_data[merge_data$eGeneID == gene, ]
      cat(sprintf("\n    Gene: %s  SNPs: %s\n", gene, nrow(gene_data)))

      if (nrow(gene_data) < 2) {
        cat("    Too few SNPs for this gene. Skipping.\n")
        failed_entries <- c(failed_entries,
                            list(data.frame(locus = locus_label, gene = gene,
                                            reason = "Too few SNPs (< 2)")))
        next
      }

      # Subset LD matrix to SNPs for this gene
      gene_snps        <- gene_data[[snp_col]]
      shared_gene_snps <- intersect(gene_snps, rownames(LD_sub))

      if (length(shared_gene_snps) < 2) {
        cat("    Insufficient LD SNP overlap for this gene. Skipping.\n")
        failed_entries <- c(failed_entries,
                            list(data.frame(locus = locus_label, gene = gene,
                                            reason = "Insufficient LD SNPs for gene")))
        next
      }

      gene_data <- gene_data[gene_data[[snp_col]] %in% shared_gene_snps, ]
      LD_gene   <- LD_sub[shared_gene_snps, shared_gene_snps, drop = FALSE]

      # Sample sizes
      GWAS_N <- 496      # Number of individuals in the GWAS
      EQTL_N <- 85      # Number of individuals in the eQTL study
  
     # n_gwas <- nrow(currloci_GWASdata)
     # n_eqtl <- nrow(currloci_eqtldata[currloci_eqtldata$eGeneID == gene, ])

      #========================
      # SECTION 12: SuSiE Fine-Mapping + Colocalization (NO FALLBACK)
      #========================

      gene_out_dir <- file.path(CurrLoci_OutDir, gene)
      dir.create(gene_out_dir, showWarnings = FALSE, recursive = TRUE)

      susie_gwas   <- NULL
      susie_eqtl   <- NULL
      susie_result <- NULL

      cat("    Running SuSiE fine-mapping on GWAS dataset...\n")
      tryCatch({
        susie_gwas <- susieR::susie_rss(
          bhat         = gene_data$beta_GWAS,
          shat         = gene_data$SE_gwas,
          R            = LD_gene,
          n            = GWAS_N,
          L            = 10,
          coverage     = 0.95,
          min_abs_corr = 0.1
        )
      }, error = function(e) {
        cat(sprintf("    SuSiE GWAS fine-mapping failed: %s\n", e$message))
        failed_entries <<- c(failed_entries,
                             list(data.frame(locus  = locus_label,
                                             gene   = gene,
                                             reason = paste("SuSiE GWAS failed:", e$message))))
      })

      if (is.null(susie_gwas)) {
        cat("    Skipping gene (SuSiE GWAS step failed).\n")
        next
      }

      cat("    Running SuSiE fine-mapping on eQTL dataset...\n")
      tryCatch({
        susie_eqtl <- susieR::susie_rss(
          bhat         = gene_data$beta_eQTL,
          shat         = gene_data$SE_eQTL,
          R            = LD_gene,
          n            = EQTL_N,
          L            = 10,
          coverage     = 0.95,
          min_abs_corr = 0.1
        )
      }, error = function(e) {
        cat(sprintf("    SuSiE eQTL fine-mapping failed: %s\n", e$message))
        failed_entries <<- c(failed_entries,
                             list(data.frame(locus  = locus_label,
                                             gene   = gene,
                                             reason = paste("SuSiE eQTL failed:", e$message))))
      })

      if (is.null(susie_eqtl)) {
        cat("    Skipping gene (SuSiE eQTL step failed).\n")
        next
      }

      cat("    Running coloc.susie colocalization...\n")
      tryCatch({
        susie_result <- coloc::coloc.susie(susie_gwas, susie_eqtl)
      }, error = function(e) {
        cat(sprintf("    coloc.susie failed: %s\n", e$message))
        failed_entries <<- c(failed_entries,
                             list(data.frame(locus  = locus_label,
                                             gene   = gene,
                                             reason = paste("coloc.susie failed:", e$message))))
      })

      if (is.null(susie_result)) {
        cat("    Skipping gene (coloc.susie step failed).\n")
        next
      }

      cat("    SuSiE colocalization completed successfully.\n")

      #========================
      # SECTION 13: Save Per-Gene SuSiE Results
      #========================

      # Save SuSiE summary table
      susie_summary <- susie_result$summary
      write.table(susie_summary,
                  file.path(gene_out_dir, "susie_coloc_summary.txt"),
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

      # Save per-SNP posterior probabilities
      if (!is.null(susie_result$results)) {
        write.table(susie_result$results,
                    file.path(gene_out_dir, "susie_coloc_results_per_snp.txt"),
                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
      }

      # Save SuSiE credible sets for GWAS and eQTL independently
      if (!is.null(susie_gwas$sets$cs)) {
        gwas_cs_df <- do.call(rbind, lapply(names(susie_gwas$sets$cs), function(cs_name) {
          data.frame(credible_set = cs_name,
                     snp_index    = susie_gwas$sets$cs[[cs_name]],
                     snp_id       = shared_gene_snps[susie_gwas$sets$cs[[cs_name]]])
        }))
        write.table(gwas_cs_df,
                    file.path(gene_out_dir, "susie_gwas_credible_sets.txt"),
                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
      }

      if (!is.null(susie_eqtl$sets$cs)) {
        eqtl_cs_df <- do.call(rbind, lapply(names(susie_eqtl$sets$cs), function(cs_name) {
          data.frame(credible_set = cs_name,
                     snp_index    = susie_eqtl$sets$cs[[cs_name]],
                     snp_id       = shared_gene_snps[susie_eqtl$sets$cs[[cs_name]]])
        }))
        write.table(eqtl_cs_df,
                    file.path(gene_out_dir, "susie_eqtl_credible_sets.txt"),
                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
      }

      # Collect colocalized signals above PP.H4 threshold
      if (!is.null(susie_summary) && "PP.H4.abf" %in% colnames(susie_summary)) {
        sig_idx <- which(susie_summary$PP.H4.abf > THR_POST_PROB)
        if (length(sig_idx) > 0) {
          cat(sprintf("    COLOCALIZATION FOUND for gene %s: %s signal(s) with PP.H4 > %.2f\n",
                      gene, length(sig_idx), THR_POST_PROB))
          sig_rows         <- susie_summary[sig_idx, ]
          sig_rows$gene    <- gene
          sig_rows$locus   <- locus_label
          sig_rows$method  <- "coloc.susie"

          if (!bool_SuSiE_Summary_DF) {
            SuSiE_Summary_DF      <- sig_rows
            bool_SuSiE_Summary_DF <- TRUE
          } else {
            SuSiE_Summary_DF <- rbind(SuSiE_Summary_DF, sig_rows)
          }
        } else {
          cat(sprintf("    No colocalization signal for gene %s (max PP.H4 = %.4f)\n",
                      gene, max(susie_summary$PP.H4.abf, na.rm = TRUE)))
        }
      }

    }  # end gene loop

  }  # end locus loop

}  # end chromosome loop

#========================
# SECTION 14: Write Final Summary Files
#========================

cat("\n\n========== Writing final output files ==========\n")

if (bool_SuSiE_Summary_DF && exists("SuSiE_Summary_DF") && nrow(SuSiE_Summary_DF) > 0) {
  write.table(SuSiE_Summary_DF, SuSiEColocSummaryFile,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  cat(sprintf("  SuSiE colocalization summary written: %s\n", SuSiEColocSummaryFile))
} else {
  cat("  No SuSiE colocalization signals found above threshold.\n")
}

# Write failed loci/gene log
if (length(failed_entries) > 0) {
  failed_df <- do.call(rbind, failed_entries)
  write.table(failed_df, SuSiEFailedLog,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  cat(sprintf("  Failed loci/genes log written: %s\n", SuSiEFailedLog))
}

cat("\n SuSiE-only pipeline complete.\n")
sink()
