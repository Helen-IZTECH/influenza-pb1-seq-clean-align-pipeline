# 02_temporal_geographic.R
# Amino-acid substitution over time and across continents at the evolvable
# contacts (757, 729). Run 01_conservation.R first (defines aa_str, year, files).
# Output: year/continent tables + Figure 2.

library(ggplot2)

# ---- 1. Continent per sequence (parsed from the file names) ----
cont_file <- toupper(sub(".*(NAMERICA|SAMERICA|EUROPE|ASIA|AFRICA|OCEANIA).*", "\\1", basename(files)))
file_n    <- sapply(files, function(f) length(fasta.seqlengths(f)))
continent <- rep(cont_file, file_n)
stopifnot(length(continent) == length(aa_str))

# ---- 2. Amino-acid frequency by year and by continent ----
for (pos in c(757, 729)) {
  a <- substr(aa_str, pos, pos)
  cat("\nPosition", pos, "- by year (%):\n")
  print(round(prop.table(table(year, a), 1) * 100, 1))
  cat("\nPosition", pos, "- by continent (%):\n")
  print(round(prop.table(table(continent, a), 1) * 100, 1))
}

# ---- 3. Figure 2: temporal substitution (757 K->R, 729 D->E) ----
yrs  <- 2009:2025
freq <- function(pos, aa) sapply(yrs, function(y)
          mean(substr(aa_str[year == as.character(y)], pos, pos) == aa) * 100)
df <- data.frame(year = rep(yrs, 2),
                 pct  = c(freq(757, "R"), freq(729, "E")),
                 site = rep(c("757 K->R", "729 D->E"), each = length(yrs)))
ggplot(df, aes(year, pct, colour = site)) +
  geom_line(linewidth = 1) + geom_point(size = 2) + ylim(0, 100) +
  labs(x = "Year", y = "% sequences with new residue", colour = NULL,
       title = "Temporal substitution at PB1-PB2 contact residues") +
  theme_classic(base_size = 13)
ggsave("fig2_temporal.png", width = 8, height = 5, dpi = 150)
