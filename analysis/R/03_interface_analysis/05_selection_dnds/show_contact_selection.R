# show_contact_selection.R
# Shows FEL / FUBAR / MEME values for the 14 (Sugiyama) and 41 (4.5A) interface sets.
# Run in the folder with fel_sites.csv, fubar_sites.csv (+ meme_sites.csv if you made it).

fel <- read.csv("fel_sites.csv")
fub <- read.csv("fubar_sites.csv")
meme_sites <- if (file.exists("meme_sites.csv")) {
  m <- read.csv("meme_sites.csv"); as.integer(m$site[as.logical(m$episodic_positive)])
} else c(361, 687, 736)   # from MEME console output

# PB1 reference amino acids (A/California/04/2009) for residue identity
REF <- paste0(
"MDVNPTLLFLKIPAQNAISTTFPYTGDPPYSHGTGTGYTMDTVNRTHQYSEKGKWTTNTETGAPQLNPID",
"GPLPEDNEPSGYAQTDCVLEAMAFLEESHPGIFENSCLETMEVVQQTRVDKLTQGRQTYDWTLNRNQPAA",
"TALANTIEVFRSNGLTANESGRLIDFLKDVMESMNKEEIEITTHFQRKRRVRDNMTKKMVTQRTIGKKKQ",
"RLNKRGYLIRALTLNTMTKDAERGKLKRRAIATPGMQIRGFVYFVETLARSICEKLEQSGLPVGGNEKKA",
"KLANVVRKMMTNSQDTEISFTITGDNTKWNENQNPRMFLAMITYITRNQPEWFRNILSMAPIMFSNKMAR",
"LGKGYMFESKRMKIRTQIPAEMLASIDLKYFNESTKKKIEKIRPLLIDGTASLSPGMMMGMFNMLSTVLG",
"VSILNLGQKKYTKTIYWWDGLQSSDDFALIVNAPNHEGIQAGVDRFYRTCKLVGINMSKKKSYINKTGTF",
"EFTSFFYRYGFVANFSMELPSFGVSGVNESADMSIGVTVIKNNMINNDLGPATAQMALQLFIKDYRYTYR",
"CHRGDTQIQTRRSFELKKLWDQTQSKVGLLVSDGGPNLYNIRNLHIPEVCLKWELMDDDYRGRLCNPLNP",
"FVSHKEIDSVNNAVVMPAHGPAKSMEYDAVATTHSWIPKRNRSILNTSQRGILEDEQMYQKCCNLFEKFF",
"PSSSYRRPVGISSMVEAMVSRARIDARVDFESGRIKKEEFSEIMKICSTIEELRRQK")
aa_of <- function(pos) substr(REF, pos, pos)
hyd <- strsplit("AVLIMFWPG","")[[1]]; pol <- strsplit("STNQCY","")[[1]]
aa_type <- function(a) ifelse(a %in% hyd,"hydrophobic", ifelse(a %in% pol,"polar","charged"))

set14 <- c(692,695,696,698,699,715,719,722,725,746,747,750,753,754)
set41 <- c(685,686,688,689,690,691,692,693,695,696,697,698,699,708,709,710,711,712,
           713,714,715,717,718,719,721,722,725,739,741,742,743,745,746,747,749,750,
           752,753,754,756,757)

show_set <- function(S, name) {
  d <- merge(fel[,c("site","alpha","beta","p_value","class")],
             fub[,c("site","purifying_post","positive_post")], by="site")
  d <- d[d$site %in% S, ]
  d$aa   <- sapply(d$site, aa_of)
  d$type <- aa_type(d$aa)
  d$FEL  <- ifelse(d$class=="negative","purifying", ifelse(d$class=="positive","POSITIVE","-"))
  d$MEME <- ifelse(d$site %in% meme_sites, "episodic+", "-")
  d <- d[order(d$site), c("site","aa","type","FEL","p_value","purifying_post","positive_post","MEME")]
  cat("\n===============  ", name, "  (", nrow(d), " residues)  ===============\n", sep="")
  print(d, row.names = FALSE)
  cat(sprintf("\nSUMMARY  FEL purifying %d | FEL positive %d | FUBAR pur>=0.9 %d | FUBAR pos>=0.9 %d | MEME episodic %d\n",
      sum(d$FEL=="purifying"), sum(d$FEL=="POSITIVE"),
      sum(d$purifying_post>=0.9), sum(d$positive_post>=0.9), sum(d$MEME=="episodic+")))
}
show_set(set14, "14 Sugiyama direct contacts")
show_set(set41, "41 (4.5A interface set)")
