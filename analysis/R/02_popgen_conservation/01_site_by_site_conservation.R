# 01_site_by_site_conservation.R
# Site-by-site nucleotide (2274 bp) and amino-acid (757 aa) conservation for the
# CLEAN H1N1 PB1 dataset. Reference: A/California/04/2009 (ACP41103.1).
# For every position: consensus, global conservation %, Shannon entropy,
# per-year conservation, year-to-year shift (|f_t+1 - f_t|), and per-continent
# conservation. Writes two styled Excel workbooks (nucleotide + amino acid).

suppressPackageStartupMessages({
  if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
  library(ape)
  library(openxlsx)
})

# ---- CONFIG -----------------------------------------------------------------
input_dir <- "data/PB1_DATASET_STANDARDIZED"   # tree containing CLEAN_ files
out_dir   <- "results"                         # where the Excel files go
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. File selection (CLEAN_ files only) ----------------------------------
all_found <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE, recursive = TRUE)
fasta_files <- all_found[grepl("^CLEAN_", basename(all_found))]
cat(sprintf("SAMPLED files excluded. Files to analyze: %d\n", length(fasta_files)))
if (length(fasta_files) == 0) stop("No FASTA files starting with CLEAN_ were found.")

# ---- 2. Reference sequence (A/California/04/2009, ACP41103.1, 757 aa) --------
REF_AA_STR <- paste0(
  "MDVNPTLLFLKIPAQNAISTTFPYTGDPPYSHGTGTGYTMDTVNRTHQYSEKGKWTTNTETGAPQLNPID",
  "GPLPEDNEPSGYAQTDCVLEAMAFLEESHPGIFENSCLETMEVVQQTRVDKLTQGRQTYDWTLNRNQPAA",
  "TALANTIEVFRSNGLTANESGRLIDFLKDVMESMNKEEIEITTHFQRKRRVRDNMTKKMVTQRTIGKKKQ",
  "RLNKRGYLIRALTLNTMTKDAERGKLKRRAIATPGMQIRGFVYFVETLARSICEKLEQSGLPVGGNEKKA",
  "KLANVVRKMMTNSQDTEISFTITGDNTKWNENQNPRMFLAMITYITRNQPEWFRNILSMAPIMFSNKMAR",
  "LGKGYMFESKRMKIRTQIPAEMLASIDLKYFNESTKKKIEKIRPLLIDGTASLSPGMMMGMFNMLSTVLG",
  "VSILNLGQKKYTKTIYWWDGLQSSDDFALIVNAPNHEGIQAGVDRFYRTCKLVGINMSKKKSYINKTGTF",
  "EFTSFFYRYGFVANFSMELPSFGVSGVNESADMSIGVTVIKNNMINNDLGPATAQMALQLFIKDYRYTYR",
  "CHRGDTQIQTRRSFELKKLWDQTQSKVGLLVSDGGPNLYNIRNLHIPEVCLKWELMDDDYRGRLCNPLNP",
  "FVSHKEIDSVNNAVVMPAHGPAKSMEYDAVATTHSWIPKRNRSILNTSQRGILEDEQMYQKCCNLFEKFF",
  "PSSSYRRPVGISSMVEAMVSRARIDARVDFESGRIKKEEFSEIMKICSTIEELRRQK"
)
REF_AA_VEC <- strsplit(REF_AA_STR, "")[[1]]

# Standard genetic code
GENCODE <- list(
  TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L",
  ATT="I", ATC="I", ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V",
  TCT="S", TCC="S", TCA="S", TCG="S", CCT="P", CCC="P", CCA="P", CCG="P",
  ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A", GCA="A", GCG="A",
  TAT="Y", TAC="Y", TAA="*", TAG="*", CAT="H", CAC="H", CAA="Q", CAG="Q",
  AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E",
  TGT="C", TGC="C", TGA="*", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R",
  AGT="S", AGC="S", AGA="R", AGG="R", GGT="G", GGC="G", GGA="G", GGG="G"
)

# ---- 3. Read sequences and translate (fast, vectorized) ---------------------
loc_map <- c("AFRICA"="AF", "ASIA"="AS", "EUROPE"="EU", "NAMERICA"="NA", "SAMERICA"="SA", "OCEANIA"="AU")
years_all <- as.character(2009:2025)
conts_all <- c("AF", "AS", "EU", "NA", "SA")

all_nt_list <- list()
all_aa_list <- list()
meta_year   <- c()
meta_cont   <- c()

cat("\nReading sequences and translating to amino acids...\n")

for (i in seq_along(fasta_files)) {
  f_path <- fasta_files[i]
  f_name <- basename(f_path)

  parts <- unlist(strsplit(tools::file_path_sans_ext(f_name), "_"))
  yr <- parts[grep("^(200[9]|201[0-9]|202[0-5])$", parts)][1]
  loc_raw <- parts[length(parts)]
  loc <- ifelse(loc_raw %in% names(loc_map), loc_map[[loc_raw]], loc_raw)

  dna <- tryCatch(read.dna(f_path, format = "fasta"), error = function(e) NULL)
  if (is.null(dna) || nrow(dna) == 0 || ncol(dna) < 2271) next

  dna_mat <- toupper(as.character(dna))[, 1:2274]
  n_s <- nrow(dna_mat)

  # Codon translation (codons 1-757)
  aa_mat <- matrix("", nrow = n_s, ncol = 757)
  for (codon_idx in 1:757) {
    c_start <- (codon_idx - 1) * 3 + 1
    codons <- paste0(dna_mat[, c_start], dna_mat[, c_start + 1], dna_mat[, c_start + 2])
    aa_mat[, codon_idx] <- unname(sapply(codons, function(c) ifelse(is.null(GENCODE[[c]]), "X", GENCODE[[c]])))
  }

  all_nt_list[[length(all_nt_list) + 1]] <- dna_mat
  all_aa_list[[length(all_aa_list) + 1]] <- aa_mat
  meta_year <- c(meta_year, rep(yr, n_s))
  meta_cont <- c(meta_cont, rep(loc, n_s))

  if (i %% 15 == 0 || i == length(fasta_files)) {
    cat(sprintf("[%d/%d] %s %s processed (N = %d)\n", i, length(fasta_files), yr, loc, n_s))
  }
}

global_nt <- do.call(rbind, all_nt_list)
global_aa <- do.call(rbind, all_aa_list)
N_TOTAL   <- nrow(global_nt)
cat(sprintf("\nTotal %d perfect PB1 sequences loaded into the matrix.\n", N_TOTAL))

# Shannon entropy over valid symbols
calc_shannon <- function(vec) {
  tb <- table(vec[vec %in% c("A","C","G","T") | (nchar(vec) == 1 & vec != "X")])
  if (length(tb) <= 1) return(0)
  p <- tb / sum(tb)
  -sum(p * log2(p))
}

# ---- 4. TABLE 1: nucleotide conservation (1-2274 bp) ------------------------
cat("\nTable 1: nucleotide conservation and yearly change (2274 bp)...\n")
nt_results <- list()

for (pos in 1:2274) {
  col_bases <- global_nt[, pos]
  col_valid <- col_bases[col_bases %in% c("A", "C", "G", "T")]

  tb_global <- sort(table(col_valid), decreasing = TRUE)
  cons_base <- names(tb_global)[1]
  cons_pct  <- round((tb_global[1] / length(col_valid)) * 100, 2)
  shannon_h <- round(calc_shannon(col_valid), 4)

  yr_pcts <- sapply(years_all, function(y) {
    b_y <- col_bases[meta_year == y]; b_y <- b_y[b_y %in% c("A", "C", "G", "T")]
    if (length(b_y) == 0) return(NA)
    round((sum(b_y == cons_base) / length(b_y)) * 100, 2)
  })

  # Year-to-year mutational shift: |f_t+1 - f_t|
  yr_diffs <- sapply(1:(length(years_all) - 1), function(idx) {
    p1 <- yr_pcts[idx]; p2 <- yr_pcts[idx + 1]
    if (is.na(p1) || is.na(p2)) return(NA)
    round(abs(p2 - p1), 2)
  })
  names(yr_diffs) <- paste0("Shift_", years_all[1:16], "_", years_all[2:17])

  cont_pcts <- sapply(conts_all, function(k) {
    b_k <- col_bases[meta_cont == k]; b_k <- b_k[b_k %in% c("A", "C", "G", "T")]
    if (length(b_k) == 0) return(NA)
    round((sum(b_k == cons_base) / length(b_k)) * 100, 2)
  })
  names(cont_pcts) <- paste0("Cont_", conts_all, "_Pct")

  nt_results[[pos]] <- c(
    Position = pos, Consensus_Base = cons_base,
    Global_Conservation_Pct = cons_pct, Shannon_Entropy = shannon_h,
    setNames(yr_pcts, paste0("Year_", years_all, "_Pct")), yr_diffs, cont_pcts
  )
}

df_nt_final <- as.data.frame(do.call(rbind, nt_results), stringsAsFactors = FALSE)
for (col in names(df_nt_final)) {
  if (col != "Consensus_Base") df_nt_final[[col]] <- as.numeric(df_nt_final[[col]])
}

# ---- 5. TABLE 2: amino-acid conservation (1-757 aa) -------------------------
cat("Table 2: amino-acid conservation and mutational shift (757 aa)...\n")
aa_results <- list()

for (pos in 1:757) {
  ref_aa  <- REF_AA_VEC[pos]
  col_aas <- global_aa[, pos]
  col_valid <- col_aas[col_aas != "" & col_aas != "X"]

  tb_global <- sort(table(col_valid), decreasing = TRUE)
  cons_aa   <- names(tb_global)[1]
  cons_pct  <- round((tb_global[1] / length(col_valid)) * 100, 2)

  ref_match_pct <- round((sum(col_valid == ref_aa) / length(col_valid)) * 100, 2)
  shannon_h     <- round(calc_shannon(col_valid), 4)

  # Major variant / escape substitution
  var_str <- "NONE (100% fixed)"
  if (length(tb_global) > 1) {
    alt_aa  <- names(tb_global)[2]
    alt_pct <- round((tb_global[2] / length(col_valid)) * 100, 2)
    var_str <- sprintf("%s -> %s (%s%%)", ref_aa, alt_aa, alt_pct)
  }

  yr_pcts <- sapply(years_all, function(y) {
    a_y <- col_aas[meta_year == y]; a_y <- a_y[a_y != "" & a_y != "X"]
    if (length(a_y) == 0) return(NA)
    round((sum(a_y == cons_aa) / length(a_y)) * 100, 2)
  })

  yr_diffs <- sapply(1:(length(years_all) - 1), function(idx) {
    p1 <- yr_pcts[idx]; p2 <- yr_pcts[idx + 1]
    if (is.na(p1) || is.na(p2)) return(NA)
    round(abs(p2 - p1), 2)
  })
  names(yr_diffs) <- paste0("Shift_", years_all[1:16], "_", years_all[2:17])

  cont_pcts <- sapply(conts_all, function(k) {
    a_k <- col_aas[meta_cont == k]; a_k <- a_k[a_k != "" & a_k != "X"]
    if (length(a_k) == 0) return(NA)
    round((sum(a_k == cons_aa) / length(a_k)) * 100, 2)
  })
  names(cont_pcts) <- paste0("Cont_", conts_all, "_Pct")

  aa_results[[pos]] <- c(
    AA_Position = pos, Ref_AA_Cal09 = ref_aa, Consensus_AA = cons_aa,
    Ref_Match_Pct = ref_match_pct, Global_Conservation_Pct = cons_pct,
    Shannon_Entropy = shannon_h, Major_Substitution = var_str,
    setNames(yr_pcts, paste0("Year_", years_all, "_Pct")), yr_diffs, cont_pcts
  )
}

df_aa_final <- as.data.frame(do.call(rbind, aa_results), stringsAsFactors = FALSE)
for (col in names(df_aa_final)) {
  if (!col %in% c("Ref_AA_Cal09", "Consensus_AA", "Major_Substitution")) {
    df_aa_final[[col]] <- as.numeric(df_aa_final[[col]])
  }
}

# ---- 6. Write styled Excel workbooks ----------------------------------------
cat("\nFormatting and saving Excel files...\n")
wb_style_header <- createStyle(fontSize = 11, fontColour = "#FFFFFF", fgFill = "#1F4E78",
                               textDecoration = "bold", halign = "center", valign = "center")
wb_style_body   <- createStyle(halign = "center", valign = "center")

wb_nt <- createWorkbook()
addWorksheet(wb_nt, "PB1_Nucleotide_Conservation")
writeData(wb_nt, "PB1_Nucleotide_Conservation", df_nt_final)
addStyle(wb_nt, "PB1_Nucleotide_Conservation", wb_style_header, rows = 1, cols = 1:ncol(df_nt_final), gridExpand = TRUE)
addStyle(wb_nt, "PB1_Nucleotide_Conservation", wb_style_body, rows = 2:(nrow(df_nt_final)+1), cols = 1:ncol(df_nt_final), gridExpand = TRUE)
setColWidths(wb_nt, "PB1_Nucleotide_Conservation", cols = 1:ncol(df_nt_final), widths = "auto")
file_nt <- file.path(out_dir, "PB1_SITE_BY_SITE_NUCLEOTIDE_CONSERVATION.xlsx")
saveWorkbook(wb_nt, file_nt, overwrite = TRUE)

wb_aa <- createWorkbook()
addWorksheet(wb_aa, "PB1_AminoAcid_Conservation")
writeData(wb_aa, "PB1_AminoAcid_Conservation", df_aa_final)
addStyle(wb_aa, "PB1_AminoAcid_Conservation", wb_style_header, rows = 1, cols = 1:ncol(df_aa_final), gridExpand = TRUE)
addStyle(wb_aa, "PB1_AminoAcid_Conservation", wb_style_body, rows = 2:(nrow(df_aa_final)+1), cols = 1:ncol(df_aa_final), gridExpand = TRUE)
setColWidths(wb_aa, "PB1_AminoAcid_Conservation", cols = 1:ncol(df_aa_final), widths = "auto")
file_aa <- file.path(out_dir, "PB1_SITE_BY_SITE_AMINOACID_CONSERVATION.xlsx")
saveWorkbook(wb_aa, file_aa, overwrite = TRUE)

cat("\n======================================================================\n")
cat("Analysis complete. Outputs:\n")
cat("1. Nucleotide table (2274 positions): ", file_nt, "\n")
cat("2. Amino-acid table (757 positions) : ", file_aa, "\n")
cat("======================================================================\n")
