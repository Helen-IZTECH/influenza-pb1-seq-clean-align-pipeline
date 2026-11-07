# combined_selection_figure.R
# Three-panel selection figure for the PB1-PB2 interface (FEL/FUBAR/MEME).
#   A) whole-protein selection landscape (MEME episodic sites ringed)
#   B) PB1 5 A structural interface set (48 residues, from binding.site on 3A1G)
#   C) Sugiyama direct-contact pocket (14 residues, subset of B)
# Inputs (from run_selection.sh): fubar_sites.csv, fel_sites.csv
# Output: FIG_selection_pb1_pb2.png / .pdf

for (p in c("ggplot2","patchwork")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(ggplot2); library(patchwork)

# setwd("C:/Users/HP/Desktop")   # <-- folder containing the CSV files

fub <- read.csv("fubar_sites.csv")
fel <- read.csv("fel_sites.csv")
fub$score <- fub$positive_post - fub$purifying_post
fub$class <- with(fub, ifelse(positive_post>=0.9,"Positive (FUBAR)",
                       ifelse(purifying_post>=0.9,"Purifying","Low signal")))
cols <- c("Purifying"="#2166AC","Positive (FUBAR)"="#D62728","Low signal"="#BDBDBD")

meme_pos  <- c(361, 687, 736)                      # MEME episodic-positive sites
invariant <- fel$site[fel$alpha==0 & fel$beta==0]  # 0 substitutions -> no selection signal
# PB1 residues within 5 A of PB2 in PDB 3A1G (binding.site, PB1-side only):
set48 <- c(685,686,688,689,690,691,692,693,694,695,696,697,698,699,700,701,708,
           709,710,711,712,713,714,715,717,718,719,721,722,725,726,729,734,735,
           739,741,742,743,745,746,747,749,750,752,753,754,756,757)
sug14 <- c(692,695,696,698,699,715,719,722,725,746,747,750,753,754)

# ---- Panel A ----
memedf <- fub[fub$site %in% meme_pos,]
pA <- ggplot(fub, aes(site, score)) +
  annotate("rect", xmin=678, xmax=757, ymin=-Inf, ymax=Inf, fill="#FDE68A", alpha=0.45) +
  geom_hline(yintercept=0, color="grey55", linewidth=0.3) +
  geom_point(aes(color=class), size=0.9, alpha=0.5) +
  geom_point(data=memedf, shape=1, size=4.2, stroke=1.1, color="black") +
  geom_text(data=memedf, aes(label=site), vjust=-1.1, size=3, fontface="bold") +
  scale_color_manual(values=cols, name=NULL) +
  scale_y_continuous(limits=c(-1.05,1.2), breaks=c(-1,-0.5,0,0.5,1)) +
  annotate("text", x=717, y=0.5, label="PB1-PB2\ninterface", size=3, color="#92400e", lineheight=0.9) +
  labs(title="A) Selection across PB1 (each dot = 1 codon)",
       subtitle="Lower = purifying, higher = positive. Red = pervasive positive (FUBAR); ring = MEME episodic.",
       x="PB1 codon position", y="selection score") +
  theme_minimal(base_size=11) +
  theme(legend.position="bottom", plot.title=element_text(face="bold"), panel.grid.minor=element_blank()) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1)))

pocket_panel <- function(resset, title, sub) {
  d <- fub[fub$site %in% resset,]
  d$pc <- with(d, ifelse(site %in% invariant,"Invariant (0 subs)",
                  ifelse(purifying_post>=0.9,"Purifying (>=0.9)","Tolerant")))
  pcols <- c("Purifying (>=0.9)"="#2166AC","Invariant (0 subs)"="#1B7837","Tolerant"="#E66101")
  ggplot(d, aes(factor(site, levels=sort(unique(site))), purifying_post, color=pc)) +
    geom_segment(aes(xend=factor(site, levels=sort(unique(site))), yend=0), linewidth=0.9) +
    geom_point(size=2.8) +
    geom_hline(yintercept=0.9, linetype="dashed", color="grey45") +
    scale_color_manual(values=pcols, name=NULL) + scale_y_continuous(limits=c(0,1.03)) +
    labs(title=title, subtitle=sub, x="PB1 residue", y="FUBAR purifying post.") +
    theme_minimal(base_size=11) +
    theme(legend.position="bottom", plot.title=element_text(face="bold"),
          panel.grid.minor=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, size=7))
}
pB <- pocket_panel(set48, "B) PB1 5A interface set (48 residues)", "Broad structural contact set. None positive.")
pC <- pocket_panel(sug14, "C) Sugiyama direct-contact pocket (14 residues)", "Literature direct contacts (subset of B). None positive.")

fig <- pA / pB / pC + plot_layout(heights=c(1,1,1))
ggsave("FIG_selection_pb1_pb2.png", fig, width=11, height=13, dpi=300)
ggsave("FIG_selection_pb1_pb2.pdf", fig, width=11, height=13)
cat("Saved FIG_selection_pb1_pb2.png / .pdf in", getwd(), "\n")
