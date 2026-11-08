library(ggplot2); library(patchwork)
d <- read.csv("ddg_vs_selection.csv")
d$class <- factor(d$class, levels=c("purifying","invariant","tolerant"))
cols <- c("purifying"="#2166AC","invariant"="#1B7837","tolerant"="#E66101")
key <- c(692,695,715,750,719,753,754,746,747)
pA <- ggplot(d, aes(class, mean_ddG, color=class)) +
  geom_boxplot(outlier.shape=NA, width=0.5) +
  geom_jitter(width=0.12, size=2.4, alpha=0.8) +
  scale_color_manual(values=cols, guide="none") +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  labs(title="A) Binding destabilization by selection class",
       subtitle="Higher = mutations more destabilizing to PB1-PB2 binding",
       x="selection class", y="mean ddG binding (EvoEF2)") +
  theme_minimal(base_size=12) + theme(plot.title=element_text(face="bold"))
pB <- ggplot(d, aes(purifying_post, mean_ddG, color=class)) +
  geom_point(size=2.6) +
  ggrepel::geom_text_repel(data=subset(d, pos %in% key), aes(label=pos), size=3, show.legend=FALSE) +
  scale_color_manual(values=cols, name=NULL) +
  labs(title="B) ddG vs evolutionary conservation",
       subtitle="Spearman rho = 0.19 (weak); tolerated sites cluster low",
       x="FUBAR purifying posterior", y="mean ddG binding") +
  theme_minimal(base_size=12) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave("FIG_ddg.png", pA+pB+plot_layout(widths=c(1,1.2)), width=12, height=5.5, dpi=160)
cat("ok\n")
