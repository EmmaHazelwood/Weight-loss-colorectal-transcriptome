library(ggplot2)
library(data.table)
library(ggrepel)
library(ggrepel)
library(ggforestplot)
library(dplyr)

dat <- fread("MR_Results/MR_combined_results.csv")
dat<-distinct(dat)
dat$Significant<-"False"
dat$Significant[which(dat$bh_p<0.05)]<-"True"
dat$LCI<-exp(dat$b-1.96*dat$se)
dat$UCI<-exp(dat$b+1.96*dat$se)
dat$OR<-exp(dat$b)

dat$Method=factor(dat$method,levels=rev(c("Inverse variance weighted","Wald ratio")))
dat$id.outcome[which(dat$id.outcome=="overall")]<-"Overall"
dat$id.outcome[which(dat$id.outcome=="colon")]<-"Colon"
dat$id.outcome[which(dat$id.outcome=="proximal")]<-"Proximal"
dat$id.outcome[which(dat$id.outcome=="distal")]<-"Distal"
dat$id.outcome[which(dat$id.outcome=="rectal")]<-"Rectal"
dat$id.outcome[which(dat$id.outcome=="left")]<-"Left"
dat$Outcome=factor(dat$id.outcome,levels=c("Overall","Colon","Proximal","Distal","Rectal","Left"))
dat$Specific<-"Colon"
dat$Specific[which(dat$type=="eQTLGen")]<-"Blood"
dat<-distinct(dat)

dat$type[which(dat$type=="eQTLGen")]<-"Blood (eQTLGen)"
dat$type[which(dat$type=="BARC")]<-"Colon (BARCUVa-Seq)"
dat$type[which(dat$type=="GTExColon")]<-"Colon (GTEx)"
dat$Instrument<-factor(dat$type,levels=rev(c("Blood (eQTLGen)","Colon (BARCUVa-Seq)","Colon (GTEx)")))

library(patchwork)
library(grid)

downlist<-c("ABHD11", "ATP5MC2", "CHMP2A", "SMAD9")  # Alphabetical order
uplist<-c("CRTC3", "EIF4ENIF1")  # Alphabetical order

down<-dat[dat$external_gene_name %in% downlist,]
down$external_gene_name <- factor(down$external_gene_name, levels = downlist)
down$Instrument <- factor(down$Instrument, levels = c("Blood (eQTLGen)","Colon (BARCUVa-Seq)","Colon (GTEx)"))

up<-dat[dat$external_gene_name %in% uplist]
up$external_gene_name <- factor(up$external_gene_name, levels = uplist)
up$Instrument <- factor(up$Instrument, levels = c("Blood (eQTLGen)","Colon (BARCUVa-Seq)","Colon (GTEx)"))

# Combine data to determine shared axis limits
# For odds ratios, we want symmetry around 1.0 on the OR scale
# So if max OR is 1.2, min should be 1/1.2 = 0.833
combined <- rbind(down, up)
or_lower <- exp(combined$b - 1.96*combined$se)
or_upper <- exp(combined$b + 1.96*combined$se)

# Find the maximum distance from 1.0
max_ratio <- max(or_upper / 1.0, 1.0 / or_lower)
xlims_or <- c(1.0 / max_ratio, 1.0 * max_ratio)

# Create individual plots for each gene
plot_list <- list()

for(i in seq_along(downlist)) {
  gene <- downlist[i]
  show_y <- (i == 1 || i == 3)  # Show y-axis for first column only (plots 1 and 3)
  
  plot_list[[gene]] <- ggforestplot::forestplot(
    df = down[down$external_gene_name == gene,],
    name = Outcome,
    estimate = b,
    pvalue = bh_p,
    psignif = 0.05,
    xlab = "",
    colour = Instrument,
    logodds = TRUE,
    xlim = xlims_or
  ) +
    scale_color_discrete(type = c("#C03830","#5B31C2","#068D94")) +
    ggtitle(gene) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = if(show_y) element_text(size = 10) else element_blank(),
      axis.text.y = if(show_y) element_text() else element_blank(),
      axis.ticks.y = if(show_y) element_line() else element_blank()
    ) +
    ylab(if(show_y) "Colorectal cancer outcome" else "")
}

for(i in seq_along(uplist)) {
  gene <- uplist[i]
  
  plot_list[[gene]] <- ggforestplot::forestplot(
    df = up[up$external_gene_name == gene,],
    name = Outcome,
    estimate = b,
    pvalue = bh_p,
    psignif = 0.05,
    xlab = "",
    colour = Instrument,
    logodds = TRUE,
    xlim = xlims_or
  ) +
    scale_color_discrete(type = c("#C03830","#5B31C2","#068D94")) +
    ggtitle(gene) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# Extract shared legend
p_legend <- ggforestplot::forestplot(
  df = down,
  name = Outcome,
  estimate = b,
  pvalue = bh_p,
  psignif = 0.05,
  colour = Instrument,
  logodds = TRUE
) +
  scale_color_discrete(type = c("#C03830","#5B31C2","#068D94")) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(title = "Instrument"))

legend <- cowplot::get_legend(p_legend)

# Create section titles
title_down <- ggplot() + 
  labs(title = "Downregulated by weight loss") +
  theme_void() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

title_up <- ggplot() + 
  labs(title = "Upregulated by weight loss") +
  theme_void() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

title_blank <- ggplot() + theme_void()

# Create x-axis label
xlab_outside <- ggplot() + 
  labs(title = "") +
  theme_void() +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

xlab_middle <- ggplot() + 
  labs(title = "Odds ratio (95% CI)") +
  theme_void() +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

# Arrange: Row 1 = titles (downreg col 1, blank col 2, upreg col 3), Rows 2-3 = plots, Row 4 = x-axis, Row 5 = legend
p <- (title_down + title_blank + title_up) /
  (plot_list[["ABHD11"]] + plot_list[["ATP5MC2"]] + plot_list[["CRTC3"]]) /
  (plot_list[["CHMP2A"]] + plot_list[["SMAD9"]] + plot_list[["EIF4ENIF1"]]) /
  (xlab_outside + xlab_middle + xlab_outside) /
  wrap_elements(legend) +
  plot_layout(heights = c(0.08, 1, 1, 0.05, 0.1), widths = c(1, 1, 1))

p
ggsave(filename="Forest Plot Combined.png", plot=last_plot(), width = 0.28*650, height = 0.28*500, units = "mm", bg = "white")