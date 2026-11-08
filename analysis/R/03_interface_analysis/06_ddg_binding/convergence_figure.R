# convergence_figure.R
# "Every interface residue is evolutionarily locked, but by different physical
#  mechanisms" -- cross-reference of three independent lines of evidence for the
#  four PB1 residues that Sugiyama et al. (2009) tested as single mutants.
#
#   Evolution  : our FEL/FUBAR selection (beta = 0 => conserved in nature)
#   Structure  : EvoEF2 binding ddG (mean over 19 subs; min = best escape sub)
#   Experiment : Sugiyama 2009 in-vivo binding + polymerase-activity assays
#
# The point: no single metric predicts conservation. V715 preserves binding
# (ddG min = -2.28; Sugiyama: binding retained) yet is the one residue whose
# mutation destroys function (-80% mRNA) -- a dynamics/signalling lock that
# binding energetics cannot see. Evolution integrates all these axes.

library(ggplot2)

axes <- c("Evolution\n(nature)", "ddG mean\n(EvoEF2)", "ddG min\n(best escape)",
          "Binding\n(Sugiyama)", "Function\n(Sugiyama)")
res  <- c("L695", "F699", "V715", "I750")

# expand.grid varies `residue` fastest, so each block of 4 is one axis column
# in residue order L695, F699, V715, I750.
df <- expand.grid(residue = factor(res, levels = res),
                  axis    = factor(axes, levels = axes),
                  stringsAsFactors = FALSE)

# label text per cell, in the same (residue-fastest, axis-slowest) order
df$label <- c(
  # L695        F699         V715         I750
  "Conserved", "Conserved", "Conserved", "Conserved",   # Evolution
  "+10.5",     "+6.2",      "+14.9",     "+13.7",        # ddG mean
  "+1.3",      "+1.6",      "-2.3",      "+2.7",         # ddG min
  "Lost",      "Weak",      "Retained",  "Lost",         # Binding
  "OK",        "Higher",    "DEAD -80%", "OK")           # Function

# lock status per cell: "lock" (green) = constrained on this axis,
#                       "free" (orange) = permissive on this axis,
#                       "escape" (red)  = a viable escape exists on this axis
df$status <- c(
  "lock","lock","lock","lock",       # Evolution: all conserved
  "lock","lock","lock","lock",       # ddG mean: all destabilising on average
  "lock","lock","escape","lock",     # ddG min: V715 has a non-destabilising sub
  "lock","lock","free","lock",       # Sugiyama binding: V715 keeps binding (not locked by binding)
  "free","free","lock","free")       # Sugiyama function: only V715 loses function

cols <- c(lock = "#1B7837", free = "#E08214", escape = "#B2182B")

p <- ggplot(df, aes(axis, residue, fill = status)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = label), color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = cols,
                    labels = c(lock = "constrained on this axis",
                               free = "permissive on this axis",
                               escape = "binding-viable escape exists"),
                    name = NULL) +
  scale_x_discrete(position = "top") +
  labs(title = "Evolution locks every PB1-PB2 interface residue — by different mechanisms",
       subtitle = "Four residues Sugiyama (2009) tested: all conserved in nature, but the physical reason differs per residue",
       x = NULL, y = NULL,
       caption = "V715: binding is retained (ddG min < 0, Sugiyama binding kept) yet function dies (-80%) — a dynamics lock binding energetics misses.") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 13),
        axis.text.x.top = element_text(size = 10, lineheight = 0.9),
        legend.position = "bottom",
        panel.grid = element_blank())

ggsave("FIG_convergence.png", p, width = 11, height = 5.2, dpi = 160)
cat("wrote FIG_convergence.png\n")
