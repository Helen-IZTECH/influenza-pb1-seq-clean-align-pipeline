# combined_selection_figure.R
# Three-method dN/dS selection figure for the PB1-PB2 interface.
# Inputs (from run_selection.sh):  fubar_sites.csv, fel_sites.csv
# Output: FIG_selection_pb1_pb2.png / .pdf
# Run in R / RStudio from the folder holding the CSV files.

for (p in c("ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(ggplot2); library(patchwork)

# setwd("C:/Users/HP/Desktop")   # <-- folder containing the CSV files

fub <- read.csv("fubar_sites.csv")
fub$is_contact <- tolower(as.character(fub$is_contact)) %in% c("true", "1")

# MEME episodic-positive sites (from meme_sites.csv / console)
meme_pos  <- c(361, 687, 736)
# contact residues with zero substitutions (FEL alpha=beta=0): report as invariant
invariant <- c(698, 722, 750, 754)

# per-site signed score: +1 positive ... -1 purifying
fub$score <- fub$positive_post - fub$purifying_post
fub$class <- with(fub, ifelse(positive_post >= 0.9, "Positive",
                       ifelse(purifying_post >= 0.9, "Purifying", "Low signal")))
cols <- c("Purifying" = "#1e3a8a", "Positive" = "#dc2626", "Low signal" = "#cbd5e1")

# ---- Panel A: whole-protein landscape ----
interface <- data.frame(xmin = 678, xmax = 757)
labs_pos  <- fub[fub$site %in% meme_pos, ]

pA <- ggplot(fub, aes(site, score)) +
  geom_rect(data = interface, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "#fde68a", alpha = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = site, yend = 0, color = class), linewidth = 0.4) +
  geom_point(data = subset(fub, is_contact), aes(site, score),
             shape = 21, color = "black", fill = "#f59e0b", size = 2.2) +
  geom_text(data = labs_pos, aes(label = site), vjust = -0.7, size = 3, fontface = "bold") +
  scale_color_manual(values = cols, name = "FUBAR class") +
  scale_y_continuous(limits = c(-1.05, 1.1)) +
  labs(title = "A) Selection landscape across PB1 (757 codons)",
       subtitle = "Down = purifying, up = positive. Yellow band = PB1-PB2 interface (678-757); orange points = 14 contact residues.",
       x = "PB1 codon position", y = "positive - purifying posterior") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ---- Panel B: contact pocket zoom ----
contact <- subset(fub, is_contact)
contact$pocket_class <- with(contact,
  ifelse(site %in% invariant, "Invariant (0 subs)",
  ifelse(purifying_post >= 0.9, "Purifying (>=0.9)", "Tolerant / low signal")))
pcols <- c("Purifying (>=0.9)" = "#1e3a8a", "Invariant (0 subs)" = "#0284c7",
           "Tolerant / low signal" = "#f59e0b")

pB <- ggplot(contact, aes(x = factor(site, levels = sort(unique(site))),
                          y = purifying_post, color = pocket_class)) +
  geom_segment(aes(xend = factor(site, levels = sort(unique(site))), yend = 0), linewidth = 1) +
  geom_point(size = 3.5) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = pcols, name = "") +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(title = "B) PB1-PB2 contact pocket (14 residues)",
       subtitle = "Purifying posterior per contact residue; none under positive selection. Dashed line = 0.9.",
       x = "PB1 contact residue", y = "FUBAR purifying posterior") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

fig <- pA / pB + plot_layout(heights = c(1, 1))
ggsave("FIG_selection_pb1_pb2.png", fig, width = 11, height = 9, dpi = 300)
ggsave("FIG_selection_pb1_pb2.pdf", fig, width = 11, height = 9)
cat("Saved FIG_selection_pb1_pb2.png / .pdf in", getwd(), "\n")
