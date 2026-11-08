# merge_ddg_selection.R
# Collapse the per-mutant ddG scan to per-residue, join with FEL/FUBAR
# selection, and classify each interface residue.
#
# Inputs  (same folder / ../05_selection_dnds output):
#   ddg_all.csv                       pos, wt, mut, ddG_bind
#   fel_sites.csv                     site, alpha, beta, p_value, ...
#   fubar_sites.csv                   site, alpha, beta, purifying_post, ...
# Output:
#   ddg_vs_selection.csv              pos, class, beta, purifying_post, mean_ddG, min_ddG
#
# Classes (PB1 numbering == alignment site):
#   invariant : FEL alpha == 0 AND beta == 0   (no substitution seen at all)
#   purifying : FUBAR purifying posterior >= 0.9  (strong purifying signal)
#   tolerant  : everything else (beta > 0, weak/no purifying signal)

ddg  <- read.csv("ddg_all.csv")
fel  <- read.csv("fel_sites.csv")
fub  <- read.csv("fubar_sites.csv")

# per-residue ddG summary
agg <- aggregate(ddG_bind ~ pos, data = ddg,
                 FUN = function(x) c(mean = mean(x), min = min(x)))
res <- data.frame(pos      = agg$pos,
                  mean_ddG = round(agg$ddG_bind[, "mean"], 3),
                  min_ddG  = round(agg$ddG_bind[, "min"], 2))

# attach selection (site == PB1 position)
res$alpha          <- fel$alpha[match(res$pos, fel$site)]
res$beta           <- fel$beta[match(res$pos, fel$site)]
res$purifying_post <- fub$purifying_post[match(res$pos, fub$site)]

res$class <- with(res, ifelse(alpha == 0 & beta == 0, "invariant",
                       ifelse(purifying_post >= 0.9, "purifying", "tolerant")))

out <- res[, c("pos", "class", "beta", "purifying_post", "mean_ddG", "min_ddG")]
out <- out[order(out$pos), ]
write.csv(out, "ddg_vs_selection.csv", row.names = FALSE)

# quick console summary
cat("n per class:\n"); print(table(out$class))
cat("\nmean ddG by class:\n")
print(round(tapply(out$mean_ddG, out$class, mean), 2))
cat(sprintf("\nSpearman(mean_ddG, purifying_post) = %.3f\n",
            cor(out$mean_ddG, out$purifying_post, method = "spearman")))
cat(sprintf("Spearman(mean_ddG, beta)           = %.3f\n",
            cor(out$mean_ddG, out$beta, method = "spearman")))
