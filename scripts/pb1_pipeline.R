## =========================================================
##  PB1 PIPELINE (STABLE)
##
##  Purpose:
##  Preprocess Influenza A PB1 nucleotide sequences for
##  downstream evolutionary and population genetic analyses.
##
##  Minimal outputs:
##   - <prefix>_final.fasta      (aligned, trimmed CDS)
##   - <prefix>_QC.csv           (QC summary table)
##
## =========================================================


## =========================================================
## (A) WORKING DIRECTORY SETUP
## =========================================================

setwd_fix <- function(path) {
  # Normalize file paths to avoid Windows backslash issues
  setwd(normalizePath(path, winslash = "/", mustWork = TRUE))
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## [EDIT-1] WORKING DIRECTORY (where input FASTA files exist)
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

setwd_fix("C:/Users/HP/Desktop/PB1_DATASET/PB1_2013/H1N1_PB1_EUROPE_2013")


## =========================================================
## (B) REQUIRED PACKAGES
## =========================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
library(Biostrings)


## =========================================================
## (C) HELPER FUNCTION: FIND FIRST ATG
## =========================================================

# Returns the position of the first ATG start codon
find_first_atg <- function(seq_char) {
  s <- toupper(gsub("U", "T", seq_char))
  n <- nchar(s)
  if (n < 3) return(NA_integer_)
  for (i in 1:(n - 2L)) {
    if (substr(s, i, i + 2L) == "ATG") return(i)
  }
  NA_integer_
}


## =========================================================
## (D) HELPER FUNCTION: TRIM FULL-GAP ALIGNMENT COLUMNS
## =========================================================

# Removes alignment columns that contain only gaps at the ends
trim_all_gap_ends <- function(aln_dnaset, gap_char = "-") {
  seqs <- as.character(aln_dnaset)
  if (length(seqs) == 0) return(aln_dnaset)
  
  m <- do.call(rbind, strsplit(seqs, split = ""))
  non_gap_counts <- colSums(m != gap_char)
  
  if (all(non_gap_counts == 0)) return(aln_dnaset)
  
  first <- min(which(non_gap_counts > 0))
  last  <- max(which(non_gap_counts > 0))
  
  m_trim  <- m[, first:last, drop = FALSE]
  seqs_tr <- apply(m_trim, 1, paste0, collapse = "")
  
  out <- DNAStringSet(seqs_tr)
  names(out) <- names(aln_dnaset)
  out
}


## =========================================================
## (E) MAIN PB1 PIPELINE
## =========================================================

run_pb1_pipeline <- function(
  input_fasta,
  output_prefix,
  mafft_path,
  mafft_threads    = 4,
  expected_cds_len = 2274,
  write_intermediate = FALSE
) {
  
  cat("===== PB1 PIPELINE START =====\n")
  cat("Input FASTA   :", input_fasta, "\n")
  cat("Output prefix :", output_prefix, "\n\n")
  
  ## -------------------------------------------------------
  ## STEP 1: READ RAW SEQUENCES AND INITIAL QC
  ## -------------------------------------------------------
  
  raw <- readDNAStringSet(input_fasta)
  raw <- DNAStringSet(toupper(raw))
  names(raw) <- gsub(" ", "_", names(raw))
  raw_char <- as.character(raw)
  
  qc <- data.frame(
    id               = names(raw),
    raw_length       = nchar(raw_char),
    invalid_char     = FALSE,
    has_gap_raw      = FALSE,
    removed_invalid  = FALSE,
    N_count          = NA_integer_,
    N_fraction       = NA_real_,
    removed_has_N    = FALSE,
    has_valid_CDS    = FALSE,
    removed_shortCDS = FALSE,
    cds_length       = NA_integer_,
    last_codon       = NA_character_,
    is_stop_last     = NA,
    aligned_length   = NA_integer_,
    gap_count        = NA_integer_,
    gap_fraction     = NA_real_,
    kept_final       = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Filter sequences containing invalid characters or gaps
  valid_pattern <- "^[ACGTURYKMSWBDHVN-]+$"
  qc$invalid_char <- !grepl(valid_pattern, raw_char)
  qc$has_gap_raw  <- grepl("-", raw_char, fixed = TRUE)
  
  removed_invalid <- qc$invalid_char | qc$has_gap_raw
  qc$removed_invalid[removed_invalid] <- TRUE
  keep_mask <- !removed_invalid
  
  cat("STEP 1: RAW\n")
  cat("  Raw sequences:", length(raw), "\n")
  cat("  Removed invalid/gap-containing:", sum(removed_invalid), "\n")
  
  
  ## -------------------------------------------------------
  ## STEP 2: REPLACE AMBIGUOUS BASES WITH N
  ## -------------------------------------------------------
  
  raw_kept   <- raw[keep_mask]
  raw_kept_N <- replaceAmbiguities(raw_kept, "N")
  raw_kept_N <- DNAStringSet(toupper(raw_kept_N))
  
  raw_kept_N_char <- as.character(raw_kept_N)
  N_count <- vapply(raw_kept_N_char, function(s) nchar(gsub("[^N]", "", s)), integer(1))
  N_frac  <- N_count / nchar(raw_kept_N_char)
  
  idx_kept <- which(keep_mask)
  qc$N_count[idx_kept]    <- N_count
  qc$N_fraction[idx_kept] <- N_frac
  
  clean <- raw_kept_N
  names(clean) <- names(raw)[keep_mask]
  
  
  ## -------------------------------------------------------
  ## STEP 3: EXTRACT FIXED-LENGTH CDS FROM FIRST ATG
  ## -------------------------------------------------------
  
  clean_char <- as.character(clean)
  n_clean    <- length(clean)
  idx_clean  <- which(keep_mask)
  
  cds_list       <- vector("list", n_clean)
  cds_keep_local <- rep(FALSE, n_clean)
  
  for (i in seq_len(n_clean)) {
    s <- clean_char[i]
    atg_pos <- find_first_atg(s)
    if (is.na(atg_pos)) {
      qc$removed_shortCDS[idx_clean[i]] <- TRUE
      next
    }
    
    end_pos <- atg_pos + expected_cds_len - 1L
    if (end_pos > nchar(s)) {
      qc$removed_shortCDS[idx_clean[i]] <- TRUE
      next
    }
    
    cds_list[[i]] <- DNAString(substr(s, atg_pos, end_pos))
    cds_keep_local[i] <- TRUE
    
    qc$has_valid_CDS[idx_clean[i]] <- TRUE
    qc$cds_length[idx_clean[i]]    <- expected_cds_len
  }
  
  cds_indices   <- which(cds_keep_local)
  cds_unaligned <- DNAStringSet(cds_list[cds_indices])
  cds_orig_idx  <- idx_clean[cds_indices]
  names(cds_unaligned) <- qc$id[cds_orig_idx]
  
  
  ## -------------------------------------------------------
  ## STEP 4: MAFFT ALIGNMENT AND TRIMMING
  ## -------------------------------------------------------
  
  cds_unaligned_fasta <- paste0(output_prefix, "_CDS_unaligned.fasta")
  writeXStringSet(cds_unaligned, cds_unaligned_fasta)
  
  cds_aligned_fasta <- paste0(output_prefix, "_CDS_aligned.fasta")
  mafft_log <- paste0(output_prefix, "_mafft.log")
  
  mafft_args <- c("--auto", "--thread", as.character(mafft_threads), cds_unaligned_fasta)
  system2(mafft_path, args = mafft_args,
          stdout = cds_aligned_fasta,
          stderr = mafft_log)
  
  cds_aligned <- readDNAStringSet(cds_aligned_fasta)
  cds_aligned_trim <- trim_all_gap_ends(cds_aligned)
  
  final_fasta <- paste0(output_prefix, "_final.fasta")
  writeXStringSet(cds_aligned_trim, final_fasta)
  
  
  ## -------------------------------------------------------
  ## STEP 5: FINAL QC AND OUTPUTS
  ## -------------------------------------------------------
  
  qc$kept_final[cds_orig_idx] <- TRUE
  qc_file <- paste0(output_prefix, "_QC.csv")
  write.csv(qc, qc_file, row.names = FALSE)
  
  cat("\n===== PB1 PIPELINE DONE =====\n")
  cat("Final alignment:", final_fasta, "\n")
  cat("QC table:", qc_file, "\n")
  
  invisible(list(
    final_alignment = cds_aligned_trim,
    cds_unaligned   = cds_unaligned,
    qc_table        = qc
  ))
}


## =========================================================
## (F) RUN PIPELINE
## =========================================================

res <- run_pb1_pipeline(
  input_fasta      = "RAW_PB1_EUROPE_2013.fasta",
  output_prefix    = "EUROPE_2013_PB1",
  mafft_path       = "C:/Users/HP/Downloads/mafft-7.526-win64-signed/mafft-win/mafft.bat",
  mafft_threads    = 4,
  expected_cds_len = 2274,
  write_intermediate = FALSE
)
