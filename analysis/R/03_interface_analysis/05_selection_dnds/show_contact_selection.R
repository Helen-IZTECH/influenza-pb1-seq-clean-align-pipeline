fel <- read.csv("fel_sites.csv")
fub <- read.csv("fubar_sites.csv")
has_meme <- file.exists("meme_sites.csv")
if (has_meme) meme <- read.csv("meme_sites.csv")

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
  d <- merge(fel[,c("site","p_value","class")], fub[,c("site","purifying_post","positive_post")], by="site")
  names(d)[names(d)=="p_value"] <- "FEL_p"
  if (has_meme) { d <- merge(d, meme[,c("site","p_value")], by="site"); names(d)[names(d)=="p_value"] <- "MEME_p" }
  d <- d[d$site %in% S, ]
  d$aa <- sapply(d$site, aa_of); d$type <- aa_type(d$aa)
  d$FEL <- ifelse(d$class=="negative","purifying", ifelse(d$class=="positive","POSITIVE","-"))
  keep <- c("site","aa","type","FEL","FEL_p","purifying_post","positive_post")
  if (has_meme) keep <- c(keep,"MEME_p")
  d <- d[order(d$site), keep]
  cat("\n===============  ", name, "  (", nrow(d), " residues)  ===============\n", sep="")
  print(d, row.names=FALSE)
  cat(sprintf("\nSUMMARY  FEL purifying %d | positive %d | FUBAR pur>=0.9 %d | pos>=0.9 %d%s\n",
      sum(d$FEL=="purifying"), sum(d$FEL=="POSITIVE"),
      sum(d$purifying_post>=0.9), sum(d$positive_post>=0.9),
      if (has_meme) sprintf(" | MEME p<=0.1: %d", sum(d$MEME_p<=0.1)) else ""))
}
show_set(set14, "14 Sugiyama")
show_set(set41, "41 (4.5A)")
