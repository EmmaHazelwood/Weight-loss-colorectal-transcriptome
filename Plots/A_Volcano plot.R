library(ggplot2)
library(data.table)
library(ggrepel)
library(ggrepel)

#Volcano plot
de<-fread("Results.csv")
de$FC<-exp(de$logFC)
de$diffexpressed <- "None"
de$diffexpressed[de$FC > 1.4 & de$P.Value<0.009369768] <- "Upregulated"
de$diffexpressed[de$FC < 0.714 & de$P.Value<0.009369768] <- "Downregulated"
de$gene<-NA
de$label<-""

# Calculate symmetrical x-axis limits
x_max <- max(abs(de$logFC), na.rm = TRUE)
x_limits <- c(-x_max, x_max)

# Filter data to points within axis limits for text labels
de_filtered <- de[abs(de$logFC) <= x_max, ]

p <- ggplot(data=de, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() + 
  theme(text = element_text(size=20), legend.position = "bottom") + 
  scale_color_manual("Differential Expression", values=c("Downregulated" = "#0072B2", "None" = "grey70", "Upregulated" = "#CC3311")) +
  xlim(x_limits)
p <- p + geom_text_repel(data=de_filtered, size=5, max.overlaps = 5000000, direction="x", nudge_x = -0.5, show.legend = FALSE)
p <- p + xlab("Natural Log Fold Change") + ylab("-log10 P Value")
p

ggsave(filename="Volcano Plot.png", plot=last_plot(),width = 1.15*240, height = 1.15*165, units = "mm",bg="white")
