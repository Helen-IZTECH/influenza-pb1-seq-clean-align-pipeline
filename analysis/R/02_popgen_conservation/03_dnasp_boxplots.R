# 03_dnasp_boxplots.R
# PB1 DnaSP metrics -> boxplots by Year and by Continent.
# Reads the workbook written by 02_dnasp_engine.R and emits one PNG per metric
# per grouping. Robust to the literal string "NA" meaning North America.

library(readxl)
library(ggplot2)

# ---- CONFIG -----------------------------------------------------------------
xlsx_path <- "results/PB1_2009_2025_DNASP_ANALYSIS.xlsx"   # from script 02
out_dir   <- "figures"
gene_name <- "PB1"
dir.create(out_dir, showWarnings = FALSE)

# ---- Read (every cell as text so "NA" is not coerced to missing) ------------
raw <- read_excel(xlsx_path, col_types = "text")

# ---- Normalize column names to safe keys, print the mapping -----------------
make_key <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)   # apostrophes, spaces, *, () -> _
  gsub("^_+|_+$", "", x)
}
old_names <- names(raw)
names(raw) <- make_key(old_names)
print(data.frame(original = old_names, key = names(raw)))  # inspect this mapping

# ---- Numeric cleaner (strips '*', takes the first numeric token) ------------
clean_num <- function(x) {
  x <- as.character(x)
  hit <- regexpr("-?[0-9]*\\.?[0-9]+", x)
  out <- rep(NA_real_, length(x))
  ok  <- !is.na(hit) & hit != -1
  if (any(ok)) out[ok] <- as.numeric(regmatches(x, hit))
  out
}

raw$year <- as.integer(clean_num(raw$year))

# Location: true blanks -> "Unknown"; literal "NA" (North America) is preserved
raw$location <- as.character(raw$location)
raw$location[is.na(raw$location) | raw$location == ""] <- "Unknown"
cont_levels <- intersect(c("AF", "AS", "EU", "NA", "SA"), unique(raw$location))
raw$location <- factor(raw$location, levels = cont_levels)

# ---- Metrics to plot: key -> axis label -------------------------------------
metrics <- c(
  syn_pol   = "Synonymous polymorphisms",
  rep_pol   = "Nonsynonymous polymorphisms",
  td        = "Tajima's D",
  fu_li_s_d = "Fu & Li's D*",
  fu_li_s_f = "Fu & Li's F*"
)
for (k in names(metrics)) raw[[k]] <- clean_num(raw[[k]])

# ---- Shared boxplot format --------------------------------------------------
make_box <- function(d, group, value, value_label, title, subtitle, file) {
  d <- d[!is.na(d[[value]]) & !is.na(d[[group]]), ]
  if (group == "year") d[[group]] <- factor(d[[group]], levels = sort(unique(d[[group]])))
  xlab <- if (group == "year") "Year" else "Continent"

  p <- ggplot(d, aes(x = .data[[group]], y = .data[[value]], color = .data[[group]])) +
    geom_boxplot(fill = "white", linewidth = 0.6, outlier.size = 1.8) +
    labs(title = title, subtitle = subtitle, x = xlab, y = value_label) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40"),
      axis.text.x   = element_text(angle = 45, hjust = 1)
    )
  ggsave(file.path(out_dir, file), p, width = 11, height = 6, dpi = 150)
}

# ---- Every metric, by year AND by continent ---------------------------------
for (k in names(metrics)) {
  lab <- metrics[[k]]
  make_box(raw, "year", k, lab,
           sprintf("%s 2009-2025 - %s", gene_name, lab), "By Year",
           sprintf("%s_%s_by_year.png", gene_name, k))
  make_box(raw, "location", k, lab,
           sprintf("%s 2009-2025 - %s", gene_name, lab), "By Continent",
           sprintf("%s_%s_by_continent.png", gene_name, k))
}

message("Done. PNGs are in: ", normalizePath(out_dir))
