# save_interface_tables_xlsx.R
# Builds styled Excel tables of FEL/FUBAR/MEME results for the 14 and 41 interface
# residues and saves PB1_interface_selection.xlsx in the working directory.
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

fel <- read.csv("fel_sites.csv")
fub <- read.csv("fubar_sites.csv")
has_meme <- file.exists("meme_sites.csv"); if (has_meme) meme <- read.csv("meme_sites.csv")

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

build <- function(S) {
  d <- merge(fel[,c("site","p_value","class")], fub[,c("site","purifying_post","positive_post")], by="site")
  names(d)[names(d)=="p_value"] <- "FEL_p"
  if (has_meme) { d <- merge(d, meme[,c("site","p_value")], by="site"); names(d)[names(d)=="p_value"] <- "MEME_p" }
  d <- d[d$site %in% S, ]
  d$aa <- sapply(d$site, aa_of); d$type <- aa_type(d$aa)
  d$FEL <- ifelse(d$class=="negative","purifying", ifelse(d$class=="positive","POSITIVE","-"))
  d$FUBAR <- ifelse(d$positive_post>=0.9,"POSITIVE", ifelse(d$purifying_post>=0.9,"purifying","-"))
  cols <- c("site","aa","type","FEL","FEL_p","FUBAR","purifying_post","positive_post")
  if (has_meme) cols <- c(cols, "MEME_p")
  d[order(d$site), cols]
}
t14 <- build(set14); t41 <- build(set41)

wb <- createWorkbook()
hdr  <- createStyle(fontColour="#FFFFFF", fgFill="#1F4E78", textDecoration="bold",
                    halign="center", valign="center", border="TopBottomLeftRight")
body <- createStyle(halign="center", border="TopBottomLeftRight")
add_sheet <- function(nm, d) {
  addWorksheet(wb, nm); writeData(wb, nm, d)
  addStyle(wb, nm, hdr, rows=1, cols=1:ncol(d), gridExpand=TRUE)
  addStyle(wb, nm, body, rows=2:(nrow(d)+1), cols=1:ncol(d), gridExpand=TRUE)
  setColWidths(wb, nm, cols=1:ncol(d), widths="auto")
  freezePane(wb, nm, firstActiveRow=2)
  pc <- which(names(d)=="purifying_post")
  conditionalFormatting(wb, nm, cols=pc, rows=2:(nrow(d)+1), type="colourScale",
                        style=c("#F8696B","#FFEB84","#63BE7B"), rule=c(0,0.5,1))
}
add_sheet("Sugiyama_14", t14)
add_sheet("Interface_41", t41)
saveWorkbook(wb, "PB1_interface_selection.xlsx", overwrite=TRUE)
cat("Saved PB1_interface_selection.xlsx in", getwd(), "\n")
