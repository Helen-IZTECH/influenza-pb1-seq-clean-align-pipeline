# 04_pp7_subtypes.R
# The PP7 binding pocket (+8 A neighbours) and whether it differs between
# PP7-sensitive and PP7-resistant influenza A subtypes.
# Needs: PDB 3A1G, and one MSA of all subtypes + reference (all_subtypes_aligned.fasta).
# Output: per-subtype consensus, discriminating positions + Figure 4.

library(bio3d)
library(Biostrings)
library(ggplot2)

# ---- 1. PP7 pocket + 8 A neighbours (chain A = PB1) ----
pdb <- read.pdb("3A1G")
pp7 <- c(690, 692, 697, 703, 708, 712, 714, 720)   # PP7 contacts (Yuan et al. 2017)
near <- binding.site(pdb,
                     a.inds = atom.select(pdb, "protein", chain = "A"),
                     b.inds = atom.select(pdb, "protein", chain = "A", resno = pp7),
                     cutoff = 8)
pocket <- sort(unique(near$resno))

# ---- 2. Per-subtype consensus, anchored to the reference ----
aln     <- readDNAStringSet("all_subtypes_aligned.fasta")   # all subtypes + reference, one MSA
mat     <- as.matrix(aln)
refi    <- grep("California", rownames(mat), ignore.case = TRUE)[1]    # reference row
nongap  <- which(mat[refi, ] != "-")                                   # reference nucleotide columns
subtype <- sub(".*(H1N1|H7N9|H9N2|H5N1|H7N7).*", "\\1", rownames(mat)) # subtype from header
stopifnot(paste(mat[refi, nongap[1:3]], collapse = "") == "ATG")       # frame check

# Consensus amino acid for one residue within a group of sequences.
consensus_aa <- function(residue, rows) {
  cols   <- nongap[c(3 * residue - 2, 3 * residue - 1, 3 * residue)]
  codons <- toupper(apply(mat[rows, cols, drop = FALSE], 1, paste, collapse = ""))
  aa <- sapply(codons, function(c) if (grepl("[^ACGT]", c)) NA else GENETIC_CODE[[c]])
  aa <- aa[!is.na(aa) & aa != "*"]
  if (!length(aa)) NA else names(sort(table(aa), decreasing = TRUE))[1]
}

subs <- c("H1N1", "H7N9", "H9N2", "H5N1", "H7N7")
cons <- sapply(subs, function(st) sapply(1:757, function(r) consensus_aa(r, which(subtype == st))))
cons <- data.frame(residue = 1:757, cons)

# ---- 3. Positions that distinguish sensitive from resistant subtypes ----
sens <- c("H1N1", "H7N9", "H9N2"); res <- c("H5N1", "H7N7")
diff <- apply(cons, 1, function(row) {
  s <- unique(na.omit(row[sens])); r <- unique(na.omit(row[res]))
  length(s) == 1 && length(r) == 1 && s != r
})
cat("Sensitive-vs-resistant differences across PB1:\n")
print(cons[diff, c("residue", sens, res)])

# ---- 4. Figure 4: pocket conservation across subtypes ----
disc <- data.frame(pos = cons$residue[diff])
ggplot() +
  annotate("rect", xmin = 678, xmax = 757, ymin = -0.12, ymax = 0.12, fill = "#4a90d9", alpha = 0.15) +
  annotate("segment", x = 1, xend = 757, y = 0, yend = 0, linewidth = 1.5, colour = "grey80") +
  geom_point(aes(pocket, 0), colour = "grey45", size = 2) +
  geom_point(aes(pp7, 0),    colour = "#2b6cb0", size = 3.5) +
  geom_point(data = disc, aes(pos, 0), colour = "#d7191c", size = 5) +
  labs(x = "PB1 residue", y = NULL,
       title = "PP7 pocket is conserved across sensitive and resistant subtypes") +
  ylim(-0.6, 0.6) + theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave("fig4_pansubtype.png", width = 11, height = 3.2, dpi = 150)
