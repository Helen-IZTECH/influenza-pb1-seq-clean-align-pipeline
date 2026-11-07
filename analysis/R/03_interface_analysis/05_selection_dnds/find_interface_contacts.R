# find_interface_contacts.R
# PB1 residues in contact with PB2 in PDB 3A1G (tetramer: 2 PB1 + 2 PB2 chains).
# Definition: PB1 residues within 4.5 A of NATIVE PB2 (resno >= 1). PB2's -2..0
# tip is a non-native construct remnant and is excluded. PB1-side residues only,
# pooled over both PB1 copies. Reproduces the 41-residue interface set.
if (!requireNamespace("bio3d", quietly = TRUE)) install.packages("bio3d")
library(bio3d)

pdb <- read.pdb("3A1G")
ca  <- pdb$atom[pdb$atom$type == "ATOM" & pdb$atom$elety == "CA", ]
maxres <- tapply(ca$resno, ca$chain, max)
pb1_chains <- names(maxres)[maxres > 600]    # PB1 = high residue numbering (~678-757)
pb2_chains <- names(maxres)[maxres <= 600]   # PB2 = low numbering (~1-37)
cat("PB1 chains:", pb1_chains, "| PB2 chains:", pb2_chains, "\n")

CUTOFF <- 4.5
get_contacts <- function(pb1_ch) {
  a  <- atom.select(pdb, chain = pb1_ch)
  b  <- atom.select(pdb, chain = pb2_chains, resno = 1:40)   # native PB2 only
  bs <- binding.site(pdb, a.inds = a, b.inds = b, cutoff = CUTOFF)
  sort(unique(bs$resno[bs$resno >= 678]))                    # PB1-side only
}
contact_residues <- sort(unique(unlist(lapply(pb1_chains, get_contacts))))
cat(length(contact_residues), "PB1 contact residues (", CUTOFF, "A):\n")
cat(contact_residues, "\n")
writeLines(paste(contact_residues, collapse = " "), "contact_residues.txt")
