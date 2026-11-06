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
# Label the two genes followed up in vitro
de$label <- ifelse(de$Symbol %in% c("ABHD11", "CHMP2A") & de$diffexpressed != "None", de$Symbol, "")

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
labelled <- de_filtered[de_filtered$label != "", ]
# Labels are placed in the empty space to the left of the downregulated genes
p <- p + geom_point(data = labelled, shape = 21, size = 3.5, stroke = 1, colour = "black", fill = NA, show.legend = FALSE) +
  geom_text_repel(data = labelled, size = 7, fontface = "italic", colour = "black",
                  nudge_x = ifelse(labelled$label == "ABHD11", -0.35, -0.7), nudge_y = ifelse(labelled$label == "ABHD11", 0.5, -1.4),
                  min.segment.length = 0, segment.colour = "black", show.legend = FALSE)
p <- p + xlab("log2 fold change") + ylab("-log10 P Value")
p

ggsave(filename="Volcano Plot.png", plot=p, width = 1.15*240, height = 1.15*135, units = "mm",bg="white")
