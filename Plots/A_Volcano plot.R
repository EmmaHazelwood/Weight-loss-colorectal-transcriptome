library(ggplot2)
library(data.table)
library(ggrepel)

#Volcano plot
de<-fread("Results.csv")
# P-value cut-off: smallest 2% of all P-values
p_cut <- sort(de$P.Value)[round(0.02*nrow(de))]
de$FC<-2^de$logFC
de$diffexpressed <- "None"
de$diffexpressed[de$logFC > log2(1.3) & de$P.Value<=p_cut] <- "Upregulated"
de$diffexpressed[de$logFC < -log2(1.3) & de$P.Value<=p_cut] <- "Downregulated"
de$gene<-NA
de$label<-""

# Symmetrical x-axis limits
x_max <- max(abs(de$logFC), na.rm = TRUE)
x_limits <- c(-x_max, x_max)

de_filtered <- de[abs(de$logFC) <= x_max, ]

p <- ggplot(data=de, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() + 
  theme(text = element_text(size=20), legend.position = "bottom") + 
  scale_color_manual("Differential Expression", values=c("Downregulated" = "#0072B2", "None" = "grey70", "Upregulated" = "#CC3311")) +
  xlim(x_limits)
p <- p + geom_text_repel(data=de_filtered, size=5, max.overlaps = 5000000, direction="x", nudge_x = -0.5, show.legend = FALSE)
p <- p + xlab("log2 fold change") + ylab("-log10 P Value")
p

ggsave(filename="Volcano Plot.png", plot=last_plot(),width = 1.15*240, height = 1.15*165, units = "mm",bg="white")
