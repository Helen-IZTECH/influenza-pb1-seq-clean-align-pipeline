# find_5A_contacts.R
# PB1 residues within 5 A of PB2 in PDB 3A1G (the structural interface set).
# binding.site returns residues on BOTH chains; keep the PB1 side (resno >= 678).
if (!requireNamespace("bio3d", quietly = TRUE)) install.packages("bio3d")
library(bio3d)

pdb <- read.pdb("3A1G")   # downloads the structure
ca <- pdb$atom[pdb$atom$type == "ATOM" & pdb$atom$elety == "CA", ]
chain_med  <- tapply(ca$resno, ca$chain, median)
pb1_chains <- names(chain_med)[chain_med > 300]
pb2_chains <- names(chain_med)[chain_med <= 300]

a.inds <- atom.select(pdb, chain = pb1_chains)
b.inds <- atom.select(pdb, chain = pb2_chains)
bs <- binding.site(pdb, a.inds = a.inds, b.inds = b.inds, cutoff = 5)

contact_pb1 <- sort(unique(bs$resno[bs$resno >= 678]))   # PB1-side contacts
cat(length(contact_pb1), "PB1 contact residues:\n"); cat(contact_pb1, "\n")
writeLines(paste(contact_pb1, collapse = " "), "contact_residues_5A.txt")
