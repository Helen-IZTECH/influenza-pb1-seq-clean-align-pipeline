# combined_selection_figure.R
# Selection figure for the PB1-PB2 interface (FEL/FUBAR/MEME).
# Inputs (from run_selection.sh): fubar_sites.csv  (+ fel_sites.csv, meme_sites.csv for reference)
# Output: FIG_selection_pb1_pb2.png / .pdf
# Run in R / RStudio from the folder holding the CSV files.

for (p in c("ggplot2","patchwork")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(ggplot2); library(patchwork)

# setwd("C:/Users/HP/Desktop")   # <-- folder containing the CSV files

fub <- read.csv("fubar_sites.csv")
fub$is_contact <- tolower(as.character(fub$is_contact)) %in% c("true","1")

meme_pos   <- c(361, 687, 736)                                   # MEME episodic-positive sites
invariant  <- c(698, 722, 750, 754)                             # contacts with 0 substitutions (FEL alpha=beta=0)
sugiyama14 <- c(692,695,696,698,699,715,719,722,725,746,747,750,753,754)

fub$score <- fub$positive_post - fub$purifying_post
fub$class <- with(fub, ifelse(positive_post>=0.9,"Positive (FUBAR)",
                       ifelse(purifying_post>=0.9,"Purifying","Low signal")))
cols <- c("Purifying"="#2166AC","Positive (FUBAR)"="#D62728","Low signal"="#BDBDBD")
memedf <- fub[fub$site %in% meme_pos,]

# ---- Panel A: whole-protein landscape ----
pA <- ggplot(fub, aes(site, score)) +
  annotate("rect", xmin=678, xmax=757, ymin=-Inf, ymax=Inf, fill="#FDE68A", alpha=0.45) +
  geom_hline(yintercept=0, color="grey55", linewidth=0.3) +
  geom_point(aes(color=class), size=0.9, alpha=0.5) +
  geom_point(data=memedf, shape=1, size=4.2, stroke=1.1, color="black") +
  geom_text(data=memedf, aes(label=site), vjust=-1.1, size=3.2, fontface="bold") +
  scale_color_manual(values=cols, name=NULL) +
  scale_y_continuous(limits=c(-1.05,1.2), breaks=c(-1,-0.5,0,0.5,1)) +
  annotate("text", x=717, y=0.55, label="PB1-PB2\ninterface", size=3.1, color="#92400e", lineheight=0.9) +
  labs(title="A) Selection across PB1 (each dot = 1 codon)",
       subtitle="Lower = purifying, higher = positive. Red = pervasive positive (FUBAR); black ring = MEME episodic positive.",
       x="PB1 codon position", y="selection score\n(-1 purifying ... +1 positive)") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom", plot.title=element_text(face="bold"), panel.grid.minor=element_blank()) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1)))

# ---- Panel B: Sugiyama direct-contact pocket (14 residues) ----
c14 <- subset(fub, site %in% sugiyama14)
c14$pocket_class <- with(c14, ifelse(site %in% invariant,"Invariant (0 subs)",
                             ifelse(purifying_post>=0.9,"Purifying (>=0.9)","Tolerant")))
pcols <- c("Purifying (>=0.9)"="#2166AC","Invariant (0 subs)"="#1B7837","Tolerant"="#E66101")
pB <- ggplot(c14, aes(factor(site, levels=sort(unique(site))), purifying_post, color=pocket_class)) +
  geom_segment(aes(xend=factor(site, levels=sort(unique(site))), yend=0), linewidth=1.1) +
  geom_point(size=4) +
  geom_hline(yintercept=0.9, linetype="dashed", color="grey45") +
  scale_color_manual(values=pcols, name=NULL) + scale_y_continuous(limits=c(0,1.03)) +
  labs(title="B) Sugiyama direct-contact pocket (14 residues)",
       subtitle="None positive. Blue = purifying, green = fully invariant, orange = tolerant rim.",
       x="PB1 contact residue", y="FUBAR purifying posterior") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom", plot.title=element_text(face="bold"), panel.grid.minor=element_blank())

fig <- pA / pB + plot_layout(heights=c(1,1))
ggsave("FIG_selection_pb1_pb2.png", fig, width=11, height=9, dpi=300)
ggsave("FIG_selection_pb1_pb2.pdf", fig, width=11, height=9)
cat("Saved FIG_selection_pb1_pb2.png / .pdf in", getwd(), "\n")
