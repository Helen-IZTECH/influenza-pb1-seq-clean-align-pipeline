# 01_conservation.R
# PB1-PB2 contact set and per-residue conservation of H1N1pdm09 (2009-2025).
# Output: invariant core vs evolvable contact residues + Figure 1.

library(Biostrings)
library(bio3d)
library(ggplot2)

# ---- 1. Load and pool all H1N1 PB1 sequences ----
data_dir <- "PB1_DATASET"                         # <- set to your root folder
files <- list.files(data_dir, pattern = "H1N1.*PB1.*\\.fasta$",
                    recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
dna_all <- do.call(c, lapply(files, readDNAStringSet))
dna_all <- DNAStringSet(gsub("[^ACGT]", "N", as.character(dna_all)))  # mask ambiguous bases
stopifnot(length(unique(width(dna_all))) == 1)                        # must be aligned (equal length)

# ---- 2. Translate; read collection year from the header ----
aa_str <- as.character(translate(dna_all, if.fuzzy.codon = "solve"))
year   <- sub(".*/", "", sapply(strsplit(names(dna_all), "\\|"), `[`, 3))  # "A/.../2025" -> "2025"

# ---- 3. Contact set: PB1 residues within 5 A of PB2 in the structure ----
pdb <- read.pdb("3A1G")                           # chain A = PB1, chain B = PB2
bs  <- binding.site(pdb,
                    a.inds = atom.select(pdb, "protein", chain = "A"),
                    b.inds = atom.select(pdb, "protein", chain = "B"),
                    cutoff = 5)
contact_residues <- sort(unique(bs$resno))

# ---- 4. Shannon entropy helper (20 standard amino acids) ----
std <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
shannon <- function(counts) {
  tot <- sum(counts); if (tot == 0) return(NA_real_)
  p <- counts / tot; p <- p[p > 0]; -sum(p * log2(p))
}

# ---- 5. Invariant core vs evolvable (bias-robust) ----
# Max entropy across well-sampled years, so large years cannot dominate the call.
yrs <- names(which(table(year) >= 100))
max_ent <- sapply(contact_residues, function(p) {
  a <- substr(aa_str, p, p)
  max(sapply(yrs, function(y) shannon(table(factor(a[year == y], levels = std)))), na.rm = TRUE)
})
core      <- contact_residues[max_ent <  0.1]
evolvable <- contact_residues[max_ent >= 0.1]
cat("Invariant core:", length(core), "residues\n")
cat("Evolvable:     ", length(evolvable), "residues ->", paste(sort(evolvable), collapse = " "), "\n")

# ---- 6. Figure 1: conservation of each contact residue ----
df <- data.frame(residue = factor(contact_residues, levels = contact_residues[order(max_ent)]),
                 max_ent = max_ent,
                 class   = ifelse(contact_residues %in% core, "invariant core", "evolvable"))
ggplot(df, aes(residue, max_ent, fill = class)) +
  geom_col() +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("invariant core" = "#999999", "evolvable" = "#D55E00")) +
  labs(x = "PB1-PB2 contact residue", y = "Maximum within-year entropy (2009-2025)",
       fill = NULL, title = "Conservation of PB1-PB2 interface contacts") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave("fig1_contact_conservation.png", width = 10, height = 5, dpi = 150)
