# 05_domain_figures.R
# Publication figures from the master workbook (04_master_workbook.R):
#   Fig 1  dumbbell + lollipop: consensus stability vs ancestral drift per domain
#   Fig 2  domain x year variability heatmap
#   Fig 3  PB1-PB2 interface conservation + continental distribution
#   Fig 4  within-domain conservation boxplot
# Also compiles a multi-page PDF report. Labels are in English.

required_pkgs <- c("readxl", "ggplot2", "dplyr", "tidyr", "patchwork", "scales", "gridExtra", "grid")
for (pkg in required_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(readxl); library(ggplot2); library(dplyr); library(tidyr)
library(patchwork); library(scales); library(gridExtra); library(grid)

# ---- CONFIG -----------------------------------------------------------------
out_dir    <- "figures"
excel_path <- "results/PB1_CONSERVATION_MASTER_WORKBOOK.xlsx"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(excel_path)) stop("Master workbook not found: ", excel_path)

df_domain <- read_excel(excel_path, sheet = "1_Structural_Domain_Summary")
df_aa_yr  <- read_excel(excel_path, sheet = "AA_1_Yearly_Conservation")
df_aa_cnt <- read_excel(excel_path, sheet = "AA_3_Continental_Distribution")

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05), color = "#1e3a8a", hjust = 0),
      plot.subtitle = element_text(size = rel(0.85), color = "#475569", margin = margin(b = 6)),
      axis.title = element_text(face = "bold", size = rel(0.88), color = "#1e293b"),
      axis.text = element_text(size = rel(0.8), color = "#334155"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f5f9", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = rel(0.8)),
      legend.text = element_text(size = rel(0.78))
    )
}

# ---- FIGURE 1: dumbbell + lollipop ------------------------------------------
df_dom_clean <- df_domain %>%
  mutate(
    Short_Domain = trimws(gsub("Domain|_", " ", Structural_Domain_Motif)),
    Variation_Rate = 100 - Mean_Conservation_Pct
  ) %>% arrange(Mean_Conservation_Pct)
df_dom_clean$Short_Domain <- factor(df_dom_clean$Short_Domain, levels = df_dom_clean$Short_Domain)

p1a <- ggplot(df_dom_clean) +
  geom_segment(aes(y = Short_Domain, yend = Short_Domain,
                   x = Cal09_Reference_Match_Pct, xend = Mean_Conservation_Pct),
               color = "#cbd5e1", linewidth = 2.2) +
  geom_point(aes(y = Short_Domain, x = Cal09_Reference_Match_Pct, color = "Cal09 reference match (%)"), size = 3.8) +
  geom_point(aes(y = Short_Domain, x = Mean_Conservation_Pct, color = "Global consensus (%)"), size = 2.4) +
  scale_color_manual(name = "Evolutionary metric:",
    values = c("Cal09 reference match (%)" = "#ef4444", "Global consensus (%)" = "#1e3a8a")) +
  scale_x_continuous(limits = c(95.5, 100.2), breaks = seq(96, 100, 1)) +
  labs(title = "A) Consensus stability vs ancestral drift",
       subtitle = "Grey line = evolutionary drift | overlapping points = fully ancestral conservation",
       x = "Rate (%)", y = "") +
  theme_pub()

p1b <- ggplot(df_dom_clean, aes(y = Short_Domain, x = Variation_Rate)) +
  geom_segment(aes(y = Short_Domain, yend = Short_Domain, x = 0, xend = Variation_Rate), color = "#0284c7", linewidth = 1.2) +
  geom_point(color = "#0284c7", size = 3.2) +
  geom_text(aes(label = sprintf("%.2f%%", Variation_Rate)), hjust = -0.3, size = 2.8, fontface = "bold", color = "#0f172a") +
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(title = "B) Variation rate (100 - conservation %)",
       subtitle = "Up to 44-fold difference across domains", x = "Variation rate (%)", y = "") +
  theme_pub() + theme(axis.text.y = element_blank())

fig1_combined <- p1a + p1b + plot_layout(ncol = 2, widths = c(1.25, 1))

# ---- FIGURE 2: domain x year variability heatmap ----------------------------
domains_def <- list(
  "N-ext (PA-binding)"   = c(1, 35),   "Fingers 1"            = c(36, 176),
  "Beta-ribbon"          = c(177, 212),"Finger tips & Pre-A"  = c(213, 264),
  "Palm (Motif A)"       = c(265, 314),"Fingers Intervening"  = c(315, 351),
  "Beta-hairpin"         = c(352, 398),"Palm (Catalytic Core)"= c(412, 498),
  "Thumb"                = c(499, 639),"Priming loop"         = c(640, 657),
  "Thumb Terminal"       = c(658, 668),"C-ext (PB2-binding)"  = c(669, 757)
)
assign_domain <- function(pos) {
  for (d in names(domains_def)) {
    rng <- domains_def[[d]]
    if (pos >= rng[1] && pos <= rng[2]) return(d)
  }
  return("Other")
}
df_aa_yr$Domain <- sapply(df_aa_yr$AA_Position, assign_domain)
year_cols <- grep("^Year_.*_Pct$", names(df_aa_yr), value = TRUE)

df_heat <- df_aa_yr %>%
  filter(Domain != "Other") %>% group_by(Domain) %>%
  summarise(across(all_of(year_cols), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("Year_"), names_to = "Year", values_to = "Conservation") %>%
  mutate(Year = gsub("Year_|_Pct", "", Year), Variation = 100 - Conservation)

p2 <- ggplot(df_heat, aes(x = Year, y = Domain, fill = Variation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Variation)), size = 2.6, color = ifelse(df_heat$Variation > 2, "black", "white")) +
  scale_fill_gradientn(name = "Variability (%):",
    colors = c("#1e3a8a", "#0284c7", "#38bdf8", "#fde047", "#ef4444"),
    values = rescale(c(0, 0.2, 0.8, 2.5, 4.5))) +
  labs(title = "H1N1 PB1 domain-wise yearly variability heatmap (2009-2025)",
       subtitle = "Dark blue = absolute conservation (100%) | light colors = mutational drift",
       x = "Year", y = "Structural domain") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

# ---- FIGURE 3: PB1-PB2 interface + continents -------------------------------
contact_residues <- c(692, 695, 696, 698, 699, 715, 719, 722, 725, 746, 747, 750, 753, 754)
cons_col_name <- grep("Global.*Cons", names(df_aa_cnt), value = TRUE)[1]
df_pb2_interface <- df_aa_cnt %>% filter(AA_Position >= 678 & AA_Position <= 757) %>%
  mutate(Is_Contact = AA_Position %in% contact_residues, Global_Cons = .data[[cons_col_name]])

p3a <- ggplot(df_pb2_interface, aes(x = AA_Position, y = Global_Cons)) +
  geom_line(color = "#cbd5e1", linetype = "dotted", linewidth = 0.8) +
  geom_point(aes(color = Is_Contact, size = Is_Contact, alpha = Is_Contact)) +
  scale_color_manual(name = "Residue type:", values = c("TRUE" = "#dc2626", "FALSE" = "#3b82f6"),
    labels = c("TRUE" = "Sugiyama contact pocket (14 AA)", "FALSE" = "Surrounding residues (66 AA)")) +
  scale_size_manual(values = c("TRUE" = 3.5, "FALSE" = 1.8), guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  annotate("text", x = 757, y = 55, label = "K757R (52.1%)", fontface = "bold", size = 3, color = "#1e3a8a", hjust = 1) +
  annotate("text", x = 729, y = 87, label = "D729E (89.3%)", fontface = "bold", size = 3, color = "#1e3a8a", hjust = 1) +
  scale_y_continuous(limits = c(50, 102), breaks = seq(50, 100, 10)) +
  labs(title = "A) PB1-PB2 interface: binding-pocket residues",
       subtitle = "Sugiyama et al. 2009 direct-contact residues highlighted in red",
       x = "PB1 amino-acid position (678-757)", y = "Global conservation (%)") + theme_pub()

sel_res <- c(757, 729, 687, 719, 739, 736, 698, 715, 754, 750, 747, 696)
df_cont_heat <- df_pb2_interface %>% filter(AA_Position %in% sel_res) %>%
  select(AA_Position, Ref_AA_Cal09, Consensus_AA, Cont_AF_Pct, Cont_AS_Pct, Cont_EU_Pct, Cont_NA_Pct, Cont_SA_Pct) %>%
  pivot_longer(cols = starts_with("Cont_"), names_to = "Continent", values_to = "Cons_Pct") %>%
  mutate(Continent = case_when(
      Continent == "Cont_AF_Pct" ~ "Africa (AF)", Continent == "Cont_AS_Pct" ~ "Asia (AS)",
      Continent == "Cont_EU_Pct" ~ "Europe (EU)", Continent == "Cont_NA_Pct" ~ "N. America (NA)",
      Continent == "Cont_SA_Pct" ~ "S. America (SA)"),
    Site_Label = paste0(Ref_AA_Cal09, AA_Position, " (", Consensus_AA, ")"))
df_cont_heat$Site_Label <- factor(df_cont_heat$Site_Label,
  levels = unique(df_cont_heat$Site_Label[order(df_cont_heat$AA_Position, decreasing = TRUE)]))

p3b <- ggplot(df_cont_heat, aes(x = Continent, y = Site_Label, fill = Cons_Pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Cons_Pct)), size = 2.8, fontface = "bold",
            color = ifelse(df_cont_heat$Cons_Pct < 70 | df_cont_heat$Cons_Pct > 96, "white", "black")) +
  scale_fill_gradientn(name = "Continental conservation (%):",
    colors = c("#fef08a", "#38bdf8", "#1e3a8a"), limits = c(25, 100)) +
  labs(title = "B) Continental distribution of critical PB1-PB2 sites",
       subtitle = "Contact-pocket (lower rows) continental homogeneity vs surface drift", x = "", y = "") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))

fig3_combined <- p3a + p3b + plot_layout(ncol = 2, widths = c(1.3, 1))

# ---- FIGURE 4: within-domain boxplot ----------------------------------------
global_aa_cons_col <- grep("Global.*Cons", names(df_aa_yr), value = TRUE)[1]
p4 <- ggplot(filter(df_aa_yr, Domain != "Other"), aes(x = .data[[global_aa_cons_col]], y = Domain)) +
  geom_boxplot(color = "#1e293b", fill = "#f8fafc", outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(color = "#dc2626", size = 1.6, alpha = 0.5, height = 0.15) +
  scale_x_continuous(limits = c(50, 102), breaks = seq(50, 100, 10)) +
  labs(title = "Within-domain conservation distribution of PB1 (N = 757 AA)",
       subtitle = "Boxes = domain body, red points = escape mutations",
       x = "Per-residue global conservation (%)", y = "") + theme_pub()

# ---- Save PNGs --------------------------------------------------------------
cat("Saving figures...\n")
ggsave(file.path(out_dir, "FIG1_domain_dumbbell_lollipop.png"), fig1_combined, width = 11, height = 6, dpi = 300)
ggsave(file.path(out_dir, "FIG2_domain_yearly_heatmap.png"), p2, width = 12, height = 6.2, dpi = 300)
ggsave(file.path(out_dir, "FIG3_pb1_pb2_interface_continents.png"), fig3_combined, width = 12.5, height = 6.2, dpi = 300)
ggsave(file.path(out_dir, "FIG4_within_domain_boxplot.png"), p4, width = 10, height = 6, dpi = 300)

# ---- Multi-page PDF report --------------------------------------------------
cat("Compiling multi-page PDF...\n")
pdf(file.path(out_dir, "H1N1_PB1_CONSERVATION_REPORT.pdf"), width = 14, height = 8.5)
df_tab_print <- df_domain %>%
  select(Structural_Domain_Motif, AA_Position_Range, Total_AA, Mean_Conservation_Pct,
         Cal09_Reference_Match_Pct, Mean_Shannon_Entropy, Most_Variable_Hotspot) %>%
  rename("Domain / Motif" = Structural_Domain_Motif, "Coordinates" = AA_Position_Range,
         "AA" = Total_AA, "Conservation %" = Mean_Conservation_Pct,
         "Cal09 match %" = Cal09_Reference_Match_Pct, "Entropy (H)" = Mean_Shannon_Entropy,
         "Most variable hotspot" = Most_Variable_Hotspot)
table_theme <- ttheme_default(
  core = list(fg_params = list(fontsize = 9), bg_params = list(fill = c("#f8fafc", "#ffffff"))),
  colhead = list(fg_params = list(fontsize = 9.5, fontface = "bold", col = "white"),
                 bg_params = list(fill = "#1e3a8a")), padding = unit(c(4, 4), "mm"))
title_grob <- textGrob(
  "Influenza A (H1N1) PB1 16-year evolutionary and structural conservation report (2009-2025)\nTarget-interface analysis prior to molecular docking",
  gp = gpar(fontsize = 14, fontface = "bold", col = "#1e3a8a"), just = "center")
subtitle_grob <- textGrob(
  "Dataset: 82 clean population FASTA | 5 continents | Models: Sugiyama et al. 2009 & Pflug et al. 2014",
  gp = gpar(fontsize = 9.5, col = "#475569"), just = "center")
tab_grob <- tableGrob(df_tab_print, rows = NULL, theme = table_theme)
grid.arrange(title_grob, subtitle_grob, tab_grob, ncol = 1, heights = c(0.12, 0.06, 0.82))
print(fig1_combined); print(p2); print(fig3_combined); print(p4)
dev.off()

cat("\n======================================================================\n")
cat("Done. Figures and PDF written to: ", normalizePath(out_dir), "\n")
cat("======================================================================\n")
