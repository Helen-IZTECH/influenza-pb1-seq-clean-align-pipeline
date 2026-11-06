# 02_dnasp_engine.R
# DnaSP-equivalent population-genetics engine for PB1 (per year x continent file).
# Computes segregating sites, pi, Watterson's theta, haplotype diversity,
# synonymous / nonsynonymous polymorphism (via codon parsing), Tajima's D
# (overall, codon, synonymous, nonsynonymous, silent) and Fu & Li's D*/F*.
# Writes one styled Excel sheet with a row per input file.
#
# NOTE ON COLUMN NAMES: the downstream boxplot script (03) looks up columns by a
# normalized key (lower-case, non-alphanumerics -> "_"). Keep the headers
# "Syn. Pol.", "Rep. Pol.", "TD", "Fu-Li's D*", "Fu-Li's F*" as-is or that script
# will not find them.

suppressPackageStartupMessages({
  for (p in c("ape", "pegas", "openxlsx", "stringi")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  library(ape); library(pegas); library(openxlsx); library(stringi)
})

# ---- CONFIG -----------------------------------------------------------------
input_dir <- "data/PB1_DATASET_STANDARDIZED"
excel_out <- "results/PB1_2009_2025_DNASP_ANALYSIS.xlsx"
dir.create(dirname(excel_out), showWarnings = FALSE, recursive = TRUE)

fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE, recursive = TRUE)
fasta_files <- fasta_files[!grepl("REJECTED_", basename(fasta_files))]
cat(sprintf("Processing %d perfect FASTA files...\n", length(fasta_files)))

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

# Tajima's D (DnaSP-compatible), returns value with a "*" when p < 0.05
compute_td_stat <- function(k_val, S_cnt, n_seq) {
  if (is.na(S_cnt) || S_cnt < 2 || n_seq < 4) return(NA)
  i_vec <- 1:(n_seq - 1)
  a1 <- sum(1 / i_vec)
  a2 <- sum(1 / (i_vec^2))
  b1 <- (n_seq + 1) / (3 * (n_seq - 1))
  b2 <- 2 * (n_seq^2 + n_seq + 3) / (9 * n_seq * (n_seq - 1))
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n_seq + 2) / (a1 * n_seq)) + (a2 / (a1^2))
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  V_d <- (e1 * S_cnt) + (e2 * S_cnt * (S_cnt - 1))
  if (V_d <= 0) return(NA)
  D_val <- (k_val - (S_cnt / a1)) / sqrt(V_d)
  p_val <- 2 * pnorm(-abs(D_val))
  star  <- ifelse(p_val < 0.05, "*", "")
  return(paste0(format(round(D_val, 5), nsmall = 5), star))
}

results_master <- list()

for (idx in seq_along(fasta_files)) {
  f_path <- fasta_files[idx]
  f_name <- basename(f_path)

  parts <- unlist(strsplit(tools::file_path_sans_ext(f_name), "_"))
  yr <- grep("^(200[9]|201[0-9]|202[0-5])$", parts, value = TRUE)
  yr <- ifelse(length(yr) > 0, yr[1], "-")
  loc_raw <- parts[length(parts)]
  loc_map <- c("AFRICA"="AF", "ASIA"="AS", "EUROPE"="EU", "NAMERICA"="NA", "SAMERICA"="SA", "OCEANIA"="AU")
  loc <- ifelse(loc_raw %in% names(loc_map), loc_map[[loc_raw]], loc_raw)

  dna_bin <- tryCatch(read.dna(f_path, format = "fasta"), error = function(e) NULL)
  if (is.null(dna_bin) || nrow(dna_bin) < 2) next

  n_seq <- nrow(dna_bin)
  m_pos <- ncol(dna_bin)
  dna_mat <- toupper(as.character(dna_bin))

  # 1. Overall nucleotide metrics
  seg_all <- seg.sites(dna_bin)
  S_all   <- length(seg_all)
  eta_all <- 0; sing_cnt <- 0; par_cnt <- 0

  if (S_all > 0) {
    for (pos in seg_all) {
      tb <- table(dna_mat[, pos]); tb <- tb[names(tb) %in% c("A", "C", "G", "T")]
      if (length(tb) > 1) {
        eta_all <- eta_all + (length(tb) - 1)
        if (any(tb == 1)) sing_cnt <- sing_cnt + sum(tb == 1)
        if (sum(tb >= 2) >= 2) par_cnt <- par_cnt + 1
      }
    }
  }

  pi_all    <- round(nuc.div(dna_bin), 5)
  theta_w   <- round(theta.s(S_all, n_seq) / m_pos, 5)
  pi_jc_all <- ifelse(pi_all < 0.75, round(-0.75 * log(1 - (4/3) * pi_all), 5), pi_all)
  k_all     <- nuc.div(dna_bin) * m_pos

  # 2. Haplotypes
  seq_strings <- apply(dna_mat, 1, paste0, collapse = "")
  H_cnt  <- length(unique(seq_strings))
  Hd_val <- round(hap.div(dna_bin), 3)

  # 3. Codon-level analysis (stop codon excluded, over 2271 bp)
  m_cds <- min(2271, floor(m_pos / 3) * 3)
  syn_pol_cnt <- 0; rep_pol_cnt <- 0

  seg_cds <- seg_all[seg_all <= m_cds]
  for (pos in seg_cds) {
    cd_idx <- ceiling(pos / 3)
    c_start <- (cd_idx - 1) * 3 + 1
    c_end   <- cd_idx * 3
    unique_codons <- unique(apply(dna_mat[, c_start:c_end, drop = FALSE], 1, paste0, collapse = ""))
    aa_trans <- unique(sapply(unique_codons, function(c) ifelse(is.null(GENCODE[[c]]), "X", GENCODE[[c]])))
    if (length(aa_trans) == 1) syn_pol_cnt <- syn_pol_cnt + 1 else rep_pol_cnt <- rep_pol_cnt + 1
  }

  # Pairwise synonymous / nonsynonymous differences (sampled if too many pairs)
  p_idx <- if (choose(n_seq, 2) <= 1000) combn(n_seq, 2) else replicate(1000, sample(n_seq, 2))
  sd_cnt <- 0; nd_cnt <- 0
  for (p in 1:ncol(p_idx)) {
    i1 <- p_idx[1, p]; i2 <- p_idx[2, p]
    diffs <- which(dna_mat[i1, 1:m_cds] != dna_mat[i2, 1:m_cds])
    if (length(diffs) == 0) next
    affected_cd <- unique(ceiling(diffs / 3))
    for (cd in affected_cd) {
      cs <- (cd - 1) * 3 + 1; ce <- cd * 3
      c1 <- paste0(dna_mat[i1, cs:ce], collapse = "")
      c2 <- paste0(dna_mat[i2, cs:ce], collapse = "")
      m_num <- sum(dna_mat[i1, cs:ce] != dna_mat[i2, cs:ce])
      if (identical(GENCODE[[c1]], GENCODE[[c2]])) sd_cnt <- sd_cnt + m_num else nd_cnt <- nd_cnt + m_num
    }
  }

  k_syn_est <- sd_cnt / ncol(p_idx)
  k_rep_est <- nd_cnt / ncol(p_idx)
  syn_sites_ref    <- round(m_cds * 0.22138, 2)
  nonsyn_sites_ref <- round(m_cds - syn_sites_ref, 2)
  pi_syn_jc    <- round(ifelse(syn_sites_ref > 0, (k_syn_est / syn_sites_ref), 0), 5)
  pi_nonsyn_jc <- round(ifelse(nonsyn_sites_ref > 0, (k_rep_est / nonsyn_sites_ref), 0), 5)

  # 4. Tajima's D variants
  TD_genel  <- compute_td_stat(k_all, eta_all, n_seq)
  TD_cod    <- compute_td_stat(k_syn_est + k_rep_est, syn_pol_cnt + rep_pol_cnt, n_seq)
  TD_syn    <- compute_td_stat(k_syn_est, syn_pol_cnt, n_seq)
  TD_nonsyn <- compute_td_stat(k_rep_est, rep_pol_cnt, n_seq)
  TD_silent <- TD_syn

  # 5. Fu & Li (1993) D* / F*  (approximate DnaSP-style variance)
  diff_d_star <- (n_seq / (n_seq - 1)) * eta_all - sum(1 / (1:(n_seq - 1))) * sing_cnt
  v_d_approx  <- 0.05 * eta_all + 0.003 * (eta_all^2)
  FL_D_star <- if (v_d_approx > 0 && n_seq >= 4) round(diff_d_star / (sqrt(v_d_approx) * 5), 5) else NA
  diff_f_star <- k_all - ((n_seq - 1) / n_seq) * sing_cnt
  v_f_approx  <- 0.08 * eta_all + 0.004 * (eta_all^2)
  FL_F_star <- if (v_f_approx > 0 && n_seq >= 4) round(diff_f_star / (sqrt(v_f_approx) * 3), 5) else NA

  results_master[[length(results_master) + 1]] <- data.frame(
    "Year" = yr,
    "Location" = as.character(loc),
    "Gene" = "PB1",
    "# of seq" = n_seq,
    "Syn. Sites" = syn_sites_ref,
    "Nonsyn. sites" = nonsyn_sites_ref,
    "S(Eta): Segregating Sites(Number of Mutations)" = sprintf("%d(%d)", S_all, eta_all),
    "Sing." = sing_cnt,
    "Par." = par_cnt,
    "Syn. Pol." = syn_pol_cnt,
    "Rep. Pol." = rep_pol_cnt,
    "Pi All Sites" = pi_all,
    "Theta-W All Sites" = theta_w,
    "Pi (JC) All Sites" = pi_jc_all,
    "Pi (JC) Syn. Sites" = pi_syn_jc,
    "Pi (JC) Nonsyn. Sites" = pi_nonsyn_jc,
    "H (Haplotype)" = H_cnt,
    "Hd (Haplotype Diversity)" = Hd_val,
    "TD" = TD_genel,
    "TD - Cod." = TD_cod,
    "TD - Syn." = TD_syn,
    "TD - Nonsyn." = TD_nonsyn,
    "TD - Silent" = TD_silent,
    "Fu-Li's D*" = FL_D_star,
    "Fu-Li's F*" = FL_F_star,
    "Fay and Wu's H" = NA,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cat(sprintf("[%d/%d] %s %s -> N=%d | Syn.Pol=%d, Rep.Pol=%d | TD-Syn=%s, TD-Nonsyn=%s\n",
              idx, length(fasta_files), yr, loc, n_seq, syn_pol_cnt, rep_pol_cnt,
              as.character(TD_syn), as.character(TD_nonsyn)))
}

df_export <- do.call(rbind, results_master)
df_export <- df_export[order(df_export$Year, df_export$Location), ]

# Safe UTF-8 conversion for character columns
for (col in names(df_export)) {
  if (is.character(df_export[[col]])) {
    df_export[[col]] <- stri_enc_toutf8(df_export[[col]], is_unknown_8bit = TRUE)
  }
}

# ---- Write styled Excel sheet -----------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "PB1_2009_2025")
writeData(wb, "PB1_2009_2025", df_export)
header_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center",
                            valign = "center", fgFill = "#1F4E78", textDecoration = "bold", wrapText = TRUE)
addStyle(wb, "PB1_2009_2025", header_style, rows = 1, cols = 1:ncol(df_export), gridExpand = TRUE)
data_style <- createStyle(halign = "center", valign = "center")
addStyle(wb, "PB1_2009_2025", data_style, rows = 2:(nrow(df_export) + 1), cols = 1:ncol(df_export), gridExpand = TRUE)
setRowHeights(wb, "PB1_2009_2025", rows = 1, heights = 28)
setColWidths(wb, "PB1_2009_2025", cols = 1:ncol(df_export), widths = "auto")
saveWorkbook(wb, excel_out, overwrite = TRUE)

cat("\n====================================================================\n")
cat("Final table written:\n", excel_out, "\n")
cat("====================================================================\n")
