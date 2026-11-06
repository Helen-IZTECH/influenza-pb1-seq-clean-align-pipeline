# 03_qc_report.R
# Per-continent quality-control and gap-evaluation report for PB1 CDS.
# Diagnoses each sequence (gaps, N content, frameshift, premature stop) and
# writes a styled two-sheet Excel workbook: a per-file decision summary and a
# per-sequence inventory.

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# ---- CONFIG -----------------------------------------------------------------
base_dir   <- "data/PB1_DATASET_STANDARDIZED"
excel_path <- file.path(base_dir, "PB1_QC_REPORT.xlsx")

if (!dir.exists(base_dir)) stop(paste("Folder not found:", base_dir))

REF_CDS_LEN <- 2274

# Scan only CLEAN and SAMPLED files (skip RAW)
all_fastas <- list.files(base_dir, pattern = "\\.(fasta|fa)$",
                         full.names = TRUE, recursive = TRUE)
target_files <- all_fastas[grepl("CLEAN_|SAMPLED", basename(all_fastas))]
cat(sprintf("Scanning %d CLEAN/SAMPLED files...\n", length(target_files)))

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

# ---- Per-sequence diagnostic ------------------------------------------------
diagnose_seq <- function(header, seq_str) {
  gap_count <- nchar(gsub("[^-]", "", seq_str))            # gap count
  clean_seq <- gsub("-", "", seq_str)                      # raw nt (gaps removed)
  nt_len    <- nchar(clean_seq)

  n_count <- nchar(gsub("[^N]", "", clean_seq))            # N count
  n_ratio <- ifelse(nt_len > 0, round((n_count / nt_len) * 100, 2), 0)

  is_in_frame <- (nt_len %% 3 == 0)                        # length multiple of 3?

  # Premature stop check (reading the gapless sequence in frame)
  premature_stop <- FALSE
  if (nt_len >= 6) {
    num_codons <- floor(nt_len / 3)
    if (num_codons > 1) {
      codons <- substring(clean_seq, seq(1, by = 3, length.out = num_codons - 1),
                          seq(3, by = 3, length.out = num_codons - 1))
      premature_stop <- any(codons %in% c("TAA", "TAG", "TGA"))
    }
  }

  # Gap-size category
  gap_category <- "ZERO_GAP"
  if (gap_count >= 1 && gap_count <= 3) {
    gap_category <- "SMALL_GAP(1-3)"
  } else if (gap_count >= 4 && gap_count <= 15) {
    gap_category <- "MEDIUM_GAP(4-15)"
  } else if (gap_count > 15) {
    gap_category <- "LARGE_GAP(>15)"
  }

  is_perfect <- (gap_count == 0 && n_count == 0 && nt_len == REF_CDS_LEN &&
                 is_in_frame && !premature_stop)

  data.frame(
    Header          = header,
    Total_Chars     = nchar(seq_str),
    Nucleotides     = nt_len,
    Gap_Count       = gap_count,
    Gap_Category    = gap_category,
    N_Count         = n_count,
    N_Ratio         = n_ratio,
    Frameshift      = !is_in_frame,
    Premature_Stop  = premature_stop,
    Perfect         = is_perfect,
    stringsAsFactors = FALSE
  )
}

summary_list  <- list()
detailed_list <- list()

for (i in seq_along(target_files)) {
  f_path <- target_files[i]
  f_name <- basename(f_path)

  # Metadata from file name
  parts <- unlist(strsplit(tools::file_path_sans_ext(f_name), "_"))
  status <- parts[1]
  yr_match <- grep("^(200[9]|201[0-9]|202[0-5])$", parts, value = TRUE)
  yr <- ifelse(length(yr_match) > 0, yr_match[1], "UNKNOWN")
  continent <- parts[length(parts)]

  parsed <- read_fasta_fast(f_path)
  n_seq  <- length(parsed$headers)
  if (n_seq == 0) next

  diag_df <- do.call(rbind, lapply(seq_len(n_seq), function(j) {
    diagnose_seq(parsed$headers[j], parsed$sequences[j])
  }))
  diag_df$File      <- f_name
  diag_df$Year      <- yr
  diag_df$Continent <- continent
  detailed_list[[i]] <- diag_df

  summary_list[[i]] <- data.frame(
    Year              = yr,
    Continent         = continent,
    File              = f_name,
    Total_Seqs        = n_seq,
    Perfect_Count     = sum(diag_df$Perfect),
    Perfect_Pct       = round((sum(diag_df$Perfect) / n_seq) * 100, 1),
    Small_Gap_1_3     = sum(diag_df$Gap_Category == "SMALL_GAP(1-3)"),
    Small_Gap_Pct     = round((sum(diag_df$Gap_Category == "SMALL_GAP(1-3)") / n_seq) * 100, 1),
    Medium_Gap_4_15   = sum(diag_df$Gap_Category == "MEDIUM_GAP(4-15)"),
    Large_Gap_Over_15 = sum(diag_df$Gap_Category == "LARGE_GAP(>15)"),
    Seqs_With_N       = sum(diag_df$N_Count > 0),
    Frameshift_Count  = sum(diag_df$Frameshift),
    Premature_Stop    = sum(diag_df$Premature_Stop),
    stringsAsFactors  = FALSE
  )
}

df_summary  <- do.call(rbind, summary_list)
df_detailed <- do.call(rbind, detailed_list)
df_summary  <- df_summary[order(df_summary$Year, df_summary$Continent), ]

# ---- Write styled Excel workbook --------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "Per_File_Decision_Summary")
writeData(wb, "Per_File_Decision_Summary", df_summary)
addWorksheet(wb, "Per_Sequence_Inventory")
writeData(wb, "Per_Sequence_Inventory", df_detailed)

hdr_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                         fgFill = "#1F4E78", textDecoration = "bold")
addStyle(wb, "Per_File_Decision_Summary", hdr_style, rows = 1, cols = 1:ncol(df_summary), gridExpand = TRUE)
addStyle(wb, "Per_Sequence_Inventory", hdr_style, rows = 1, cols = 1:ncol(df_detailed), gridExpand = TRUE)
setColWidths(wb, "Per_File_Decision_Summary", cols = 1:ncol(df_summary), widths = "auto")
setColWidths(wb, "Per_Sequence_Inventory", cols = 1:ncol(df_detailed), widths = "auto")

saveWorkbook(wb, excel_path, overwrite = TRUE)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("QC report written:\n")
cat(excel_path, "\n")
