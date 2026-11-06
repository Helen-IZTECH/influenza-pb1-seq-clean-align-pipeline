# 04_master_workbook.R
# Builds the master conservation workbook (9 modular sheets) from the two
# site-by-site tables written by 01_site_by_site_conservation.R.
# Includes a dedicated Sugiyama et al. 2009 PB1-PB2 contact-pocket sheet.
# Column detection is pattern-based, so it tolerates small header changes.

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
  if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
  library(openxlsx); library(readxl)
})

# ---- CONFIG -----------------------------------------------------------------
results_dir     <- "results"
raw_aa_path     <- file.path(results_dir, "PB1_SITE_BY_SITE_AMINOACID_CONSERVATION.xlsx")
raw_nt_path     <- file.path(results_dir, "PB1_SITE_BY_SITE_NUCLEOTIDE_CONSERVATION.xlsx")
master_wb_path  <- file.path(results_dir, "PB1_CONSERVATION_MASTER_WORKBOOK.xlsx")

if (!file.exists(raw_aa_path)) stop("Amino-acid table not found: ", raw_aa_path)
if (!file.exists(raw_nt_path)) stop("Nucleotide table not found: ", raw_nt_path)

cat("Reading site-by-site tables...\n")
df_aa <- read_excel(raw_aa_path)
df_nt <- read_excel(raw_nt_path)

# ---- 1. Smart column detection (amino acid) ---------------------------------
cols_aa <- names(df_aa)
aa_pos_col  <- cols_aa[grep("pos", cols_aa, ignore.case = TRUE)][1]
aa_ref_aa   <- cols_aa[grep("ref_aa|ref.*cal", cols_aa, ignore.case = TRUE)][1]
aa_cons_aa  <- cols_aa[grep("consensus.*aa|cons.*aa", cols_aa, ignore.case = TRUE)][1]
aa_cons_col <- cols_aa[grep("global.*cons|cons.*pct|cons.*%", cols_aa, ignore.case = TRUE)][1]
aa_ref_col  <- cols_aa[grep("ref.*match|ref.*pct|ref.*%", cols_aa, ignore.case = TRUE)][1]
aa_ent_col  <- cols_aa[grep("entropy|shannon", cols_aa, ignore.case = TRUE)][1]
aa_mut_col  <- cols_aa[grep("major|subst|mutat", cols_aa, ignore.case = TRUE)][1]

df_aa$POS  <- as.numeric(df_aa[[aa_pos_col]])
df_aa$CONS <- as.numeric(gsub("[%,]", ".", as.character(df_aa[[aa_cons_col]])))
df_aa$REF  <- if (!is.na(aa_ref_col)) as.numeric(gsub("[%,]", ".", as.character(df_aa[[aa_ref_col]]))) else df_aa$CONS
df_aa$ENT  <- if (!is.na(aa_ent_col)) as.numeric(df_aa[[aa_ent_col]]) else NA
df_aa$MUT  <- if (!is.na(aa_mut_col)) as.character(df_aa[[aa_mut_col]]) else "NONE"

aa_core_cols  <- c(aa_pos_col, aa_ref_aa, aa_cons_aa, aa_cons_col, aa_ref_col, aa_ent_col, aa_mut_col)
aa_core_cols  <- aa_core_cols[!is.na(aa_core_cols)]
aa_year_cols  <- cols_aa[grep("^Year_.*_Pct$|^[0-9]{4}_AA_%$", cols_aa)]
aa_shift_cols <- cols_aa[grep("^Shift_", cols_aa)]
aa_cont_cols  <- cols_aa[grep("^Cont_.*_Pct$|^[A-Z]{2}_%$", cols_aa)]

# ---- 2. Smart column detection (nucleotide) ---------------------------------
cols_nt <- names(df_nt)
nt_pos_col  <- cols_nt[grep("pos", cols_nt, ignore.case = TRUE)][1]
nt_cons_col <- cols_nt[grep("consensus.*base|cons.*base|cons.*nt", cols_nt, ignore.case = TRUE)][1]
nt_glob_col <- cols_nt[grep("global.*cons|cons.*pct|cons.*%", cols_nt, ignore.case = TRUE)][1]
nt_ent_col  <- cols_nt[grep("entropy|shannon", cols_nt, ignore.case = TRUE)][1]

nt_core_cols  <- c(nt_pos_col, nt_cons_col, nt_glob_col, nt_ent_col)
nt_core_cols  <- nt_core_cols[!is.na(nt_core_cols)]
nt_year_cols  <- cols_nt[grep("^Year_.*_Pct$|^[0-9]{4}_%$", cols_nt)]
nt_shift_cols <- cols_nt[grep("^Shift_", cols_nt)]
nt_cont_cols  <- cols_nt[grep("^Cont_.*_Pct$|^[A-Z]{2}_%$", cols_nt)]

# ---- 3. SHEET 1: structural domain summary (Pflug et al. 2014 model) --------
domains_def <- list(
  "N_ext_PA_Binding"       = c(1, 35),
  "Fingers_Domain_1"       = c(36, 176),
  "Beta_ribbon"            = c(177, 212),
  "Finger_tips_and_Pre_A"  = c(213, 264),
  "Palm_Domain_Motif_A"    = c(265, 314),
  "Fingers_Intervening"    = c(315, 351),
  "Beta_hairpin"           = c(352, 398),
  "Palm_Catalytic_Core"    = c(412, 498),  # Motifs B, C(SDD), D, E
  "Thumb_Domain"           = c(499, 639),
  "Priming_loop"           = c(640, 657),
  "Thumb_Terminal"         = c(658, 668),
  "C_ext_PB2_Binding"      = c(669, 757)   # PB2-binding region
)

summary_rows <- list()
for (d_name in names(domains_def)) {
  rng <- domains_def[[d_name]]
  sub_df <- df_aa[df_aa$POS >= rng[1] & df_aa$POS <= rng[2] & !is.na(df_aa$CONS), ]
  if (nrow(sub_df) == 0) next
  mean_cons <- round(mean(sub_df$CONS, na.rm = TRUE), 2)
  mean_ref  <- round(mean(sub_df$REF, na.rm = TRUE), 2)
  mean_ent  <- round(mean(sub_df$ENT, na.rm = TRUE), 4)
  min_idx <- which.min(sub_df$CONS)
  hotspot_info <- if (length(min_idx) > 0 && !is.na(min_idx)) {
    sprintf("AA %d (%s, %.1f%%)", sub_df$POS[min_idx], sub_df$MUT[min_idx], sub_df$CONS[min_idx])
  } else "NONE (100% fixed)"

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Structural_Domain_Motif   = gsub("_", " ", d_name),
    AA_Position_Range         = sprintf("%d - %d", rng[1], rng[2]),
    Total_AA                  = nrow(sub_df),
    Mean_Conservation_Pct     = mean_cons,
    Cal09_Reference_Match_Pct = mean_ref,
    Mean_Shannon_Entropy      = mean_ent,
    Most_Variable_Hotspot     = hotspot_info,
    stringsAsFactors = FALSE
  )
}
table1_domains <- do.call(rbind, summary_rows)

# ---- 4. SHEET 2: variable sites (< 98.5% conserved) -------------------------
table2_variable <- df_aa[df_aa$CONS < 98.5 & !is.na(df_aa$CONS), ]
table2_variable <- table2_variable[order(table2_variable$CONS), ]
sel_sample_yrs <- aa_year_cols[grep("2009|2015|2020|2025", aa_year_cols)]
cols_tab2 <- c(aa_core_cols, sel_sample_yrs, aa_shift_cols[1:min(2, length(aa_shift_cols))])
table2_compact <- table2_variable[, intersect(cols_tab2, names(table2_variable))]

# ---- 5. SHEET 3: Sugiyama et al. 2009 PB1-PB2 contact pocket ----------------
cat("Building the Sugiyama et al. 2009 PB1-PB2 contact sheet...\n")
sugiyama_contact_res <- c(692, 695, 696, 698, 699, 715, 719, 722, 725, 746, 747, 750, 753, 754)

df_interface <- df_aa[df_aa$POS >= 678 & df_aa$POS <= 757, ]
df_interface$Interaction_Role <- ifelse(
  df_interface$POS %in% sugiyama_contact_res,
  "Sugiyama direct contact pocket (14 AA)",
  "Surrounding / surface residue (66 AA)"
)
# Put the 14 contact residues first, the rest in positional order
df_sugiyama_tab <- df_interface[order(!(df_interface$POS %in% sugiyama_contact_res), df_interface$POS), ]

sel_sugi_cols <- c("Interaction_Role", aa_pos_col, aa_ref_aa, aa_cons_aa,
                   aa_cons_col, aa_ref_col, aa_ent_col, aa_cont_cols, sel_sample_yrs)
sel_sugi_cols <- intersect(sel_sugi_cols, names(df_sugiyama_tab))
table3_sugiyama <- df_sugiyama_tab[, sel_sugi_cols]

# ---- 6. Remaining modular tables --------------------------------------------
tab_aa_years <- df_aa[, c(aa_core_cols, aa_year_cols)]
tab_aa_shift <- df_aa[, c(aa_core_cols[1:3], aa_shift_cols)]
tab_aa_conts <- df_aa[, c(aa_core_cols[1:3], aa_core_cols[4], aa_cont_cols)]
tab_nt_years <- df_nt[, c(nt_core_cols, nt_year_cols)]
tab_nt_shift <- df_nt[, c(nt_core_cols[1:2], nt_shift_cols)]
tab_nt_conts <- df_nt[, c(nt_core_cols[1:3], nt_cont_cols)]

# ---- 7. Build and format the workbook ---------------------------------------
cat("\nFormatting and saving the master workbook...\n")
wb <- createWorkbook()
style_hdr <- createStyle(fontSize = 11, fontColour = "#FFFFFF", fgFill = "#1F4E78",
                         textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE)
style_data   <- createStyle(halign = "center", valign = "center", fontSize = 10)
style_pocket <- createStyle(fgFill = "#FCE4D6", fontColour = "#C00000", textDecoration = "bold",
                            halign = "center", valign = "center", fontSize = 10)

add_formatted_sheet <- function(wb, sheet_name, data_df, is_summary = FALSE, is_sugiyama = FALSE) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data_df)
  addStyle(wb, sheet_name, style_hdr, rows = 1, cols = 1:ncol(data_df), gridExpand = TRUE)
  addStyle(wb, sheet_name, style_data, rows = 2:(nrow(data_df) + 1), cols = 1:ncol(data_df), gridExpand = TRUE)
  if (is_sugiyama) {
    addStyle(wb, sheet_name, style_pocket, rows = 2:15, cols = 1:ncol(data_df), gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = 1, widths = 38)
    setColWidths(wb, sheet_name, cols = 2:ncol(data_df), widths = 16)
  } else {
    setColWidths(wb, sheet_name, cols = 1:ncol(data_df), widths = "auto")
  }
  setRowHeights(wb, sheet_name, rows = 1, heights = if (is_summary) 28 else 24)
}

# Analytical summary sheets
add_formatted_sheet(wb, "1_Structural_Domain_Summary", table1_domains, is_summary = TRUE)
add_formatted_sheet(wb, "2_Variable_Sites", table2_compact)
add_formatted_sheet(wb, "3_Sugiyama_PB1_PB2_Pocket", table3_sugiyama, is_summary = TRUE, is_sugiyama = TRUE)
# Amino-acid detail sheets
add_formatted_sheet(wb, "AA_1_Yearly_Conservation", tab_aa_years)
add_formatted_sheet(wb, "AA_2_YearToYear_Shift", tab_aa_shift)
add_formatted_sheet(wb, "AA_3_Continental_Distribution", tab_aa_conts)
# Nucleotide detail sheets
add_formatted_sheet(wb, "NT_1_Yearly_Conservation", tab_nt_years)
add_formatted_sheet(wb, "NT_2_YearToYear_Shift", tab_nt_shift)
add_formatted_sheet(wb, "NT_3_Continental_Distribution", tab_nt_conts)

saveWorkbook(wb, master_wb_path, overwrite = TRUE)

cat("\n======================================================================\n")
cat("Done. 9-sheet master workbook written:\n")
cat(master_wb_path, "\n")
cat("======================================================================\n")
