# 03_linkage_phylogeny.R
# Are the two evolving residues (757R, 729E) linked? Co-occurrence, a
# year-matched test, and a phylogenetic tree coloured by residue state.
# Run 01_conservation.R first (defines aa_str, year).
# Needs a 500-sequence subsample FASTA: PB1_subsample.fasta.

library(Biostrings)
library(ape)
library(phangorn)

r757 <- substr(aa_str, 757, 757) == "R"
e729 <- substr(aa_str, 729, 729) == "E"

# ---- 1. Co-occurrence of the two variants ----
cat("Counts:\n"); print(addmargins(table(R757 = r757, E729 = e729)))

# ---- 2. Year-matched test (2024-2025, where both variants circulate) ----
# Restricting to these years controls for the temporal sweep.
yr <- year %in% c("2024", "2025")
cat("\n2024-25  729E | on 757R:", round(mean(e729[yr &  r757]) * 100, 1), "%",
    "|  729E | on 757K:",        round(mean(e729[yr & !r757]) * 100, 1), "%\n")
print(fisher.test(table(R757 = r757[yr], E729 = e729[yr])))

# ---- 3. NJ tree of a subsample; tips coloured by residue state ----
sub  <- readDNAStringSet("PB1_subsample.fasta")
tree <- NJ(dist.ml(phyDat(as.matrix(sub), type = "DNA")))
aa   <- as.character(translate(DNAStringSet(gsub("[^ACGT]", "N", as.character(sub))),
                               if.fuzzy.codon = "solve"))
names(aa) <- names(sub)
a757 <- substr(aa[tree$tip.label], 757, 757)
a729 <- substr(aa[tree$tip.label], 729, 729)

png("fig3_tree.png", width = 1400, height = 700, res = 130)
par(mfrow = c(1, 2))
plot(tree, type = "fan", show.tip.label = FALSE, main = "757 (R = orange)")
tiplabels(pch = 19, cex = 0.5, col = ifelse(a757 == "R", "#D55E00", "grey75"))
plot(tree, type = "fan", show.tip.label = FALSE, main = "729 (E = blue)")
tiplabels(pch = 19, cex = 0.5, col = ifelse(a729 == "E", "#0072B2", "grey75"))
dev.off()
