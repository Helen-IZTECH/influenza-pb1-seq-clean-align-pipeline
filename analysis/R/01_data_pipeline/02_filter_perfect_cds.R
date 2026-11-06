# 02_filter_perfect_cds.R
# Strict quality filter for PB1 coding sequences.
# Keeps only "perfect" CDS: exactly 2274 nt, ATG start, valid stop codon,
# no premature stop, no N / non-ACGT characters. Everything else is written
# to a REJECTED_ file next to the input, plus a CSV filtering report.

# ---- CONFIG -----------------------------------------------------------------
base_dir   <- "data/PB1_DATASET_STANDARDIZED"  # input tree (from script 01)
output_dir <- "data/PB1_DATASET_STANDARDIZED"  # output tree (in place is fine)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

REF_CDS_LEN <- 2274

# ---- Fast FASTA reader ------------------------------------------------------
read_fasta_fast <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  header_idx <- which(grepl("^>", lines))
  if (length(header_idx) == 0) return(list())

  headers <- sub("^>", "", lines[header_idx])
  seq_starts <- header_idx + 1
  seq_ends   <- c(header_idx[-1] - 1, length(lines))

  seqs <- character(length(header_idx))
  for (i in seq_along(header_idx)) {
    if (seq_starts[i] <= seq_ends[i]) {
      chunk <- lines[seq_starts[i]:seq_ends[i]]
      seqs[i] <- toupper(gsub("[ \t\r\n]", "", paste(chunk, collapse = "")))
    }
  }
  list(headers = headers, sequences = seqs)
}

# ---- Perfection test --------------------------------------------------------
is_sequence_perfect <- function(seq_str) {
  clean_seq <- gsub("[^A-Z]", "", seq_str)                 # drop gaps/whitespace
  if (nchar(clean_seq) != REF_CDS_LEN) return(FALSE)       # exact length
  if (grepl("[^ATCG]", clean_seq)) return(FALSE)           # no N / ambiguous
  if (substr(clean_seq, 1, 3) != "ATG") return(FALSE)      # start codon
  last_codon <- substr(clean_seq, 2272, 2274)              # stop codon
  if (!last_codon %in% c("TAA", "TAG", "TGA")) return(FALSE)
  # no premature stop among the first 757 codons
  codons <- substring(clean_seq, seq(1, by = 3, length.out = 757),
                      seq(3, by = 3, length.out = 757))
  if (any(codons %in% c("TAA", "TAG", "TGA"))) return(FALSE)
  return(TRUE)
}

# ---- FASTA writer -----------------------------------------------------------
write_fasta_clean <- function(headers, seqs, filepath) {
  if (length(headers) == 0) return()
  lines <- vector("character", length(headers) * 2)
  lines[seq(1, by = 2, length.out = length(headers))] <- paste0(">", headers)
  lines[seq(2, by = 2, length.out = length(headers))] <- seqs
  writeLines(lines, filepath)
}

# ---- Scan CLEAN_ and SAMPLED files ------------------------------------------
all_fastas <- list.files(base_dir, pattern = "\\.(fasta|fa)$",
                         full.names = TRUE, recursive = TRUE)
target_files <- all_fastas[grepl("CLEAN_|SAMPLED", basename(all_fastas))]

cat(sprintf("Filtering %d files...\n", length(target_files)))
cat(paste(rep("=", 80), collapse = ""), "\n\n")

total_kept <- 0
total_rejected <- 0
log_report <- list()

for (i in seq_along(target_files)) {
  f_path <- target_files[i]
  f_name <- basename(f_path)

  parsed <- read_fasta_fast(f_path)
  n_total <- length(parsed$headers)
  if (n_total == 0) next

  perfect_idx <- sapply(parsed$sequences, is_sequence_perfect)

  kept_headers <- parsed$headers[perfect_idx]
  kept_seqs    <- gsub("[^A-Z]", "", parsed$sequences[perfect_idx])  # gapless nt
  rej_headers  <- parsed$headers[!perfect_idx]
  rej_seqs     <- parsed$sequences[!perfect_idx]

  # Preserve the hierarchical folder structure in the output tree
  rel_dir <- dirname(sub(paste0("^", normalizePath(base_dir, winslash = "/"), "/?"), "",
                         normalizePath(f_path, winslash = "/")))
  dest_dir <- file.path(output_dir, rel_dir)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

  if (length(kept_headers) > 0) {
    write_fasta_clean(kept_headers, kept_seqs, file.path(dest_dir, f_name))
  }
  if (length(rej_headers) > 0) {
    write_fasta_clean(rej_headers, rej_seqs, file.path(dest_dir, paste0("REJECTED_", f_name)))
  }

  n_kept <- length(kept_headers)
  n_rej  <- length(rej_headers)
  total_kept <- total_kept + n_kept
  total_rejected <- total_rejected + n_rej

  log_report[[i]] <- data.frame(
    File          = f_name,
    Original      = n_total,
    Kept_Perfect  = n_kept,
    Rejected      = n_rej,
    Kept_Pct      = round((n_kept / n_total) * 100, 1),
    stringsAsFactors = FALSE
  )

  if (i %% 10 == 0 || i == length(target_files)) {
    cat(sprintf("[%d/%d] %-38s -> kept: %d | rejected: %d\n",
                i, length(target_files), f_name, n_kept, n_rej))
  }
}

df_log <- do.call(rbind, log_report)
write.csv(df_log, file.path(output_dir, "FILTERING_REPORT.csv"), row.names = FALSE)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Filtering complete.\n")
cat(sprintf("Total kept (100%% perfect): %d sequences\n", total_kept))
cat(sprintf("Total rejected: %d sequences\n", total_rejected))
cat(sprintf("Clean dataset location:\n%s\n", output_dir))
