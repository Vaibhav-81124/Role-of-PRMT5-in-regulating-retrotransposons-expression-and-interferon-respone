library(ggplot2)
library(RColorBrewer)
library(ggrepel)

df <- read.table("diff-exp-fdr.txt")

## Create a new column named diff.exp to catergorize the genes by direction of change
df$diff.exp <- "NO"
## set criteria for differential expression and replace NO with UP orr Down
df$diff.exp[df$log2FoldChange > 0.6 & df$pvalue < 0.05] <- "UP"
df$diff.exp[df$log2FoldChange < -0.6 & df$pvalue < 0.05] <- "DOWN"

## get colors
col=brewer.pal(n=8, name="Dark2")

col
#[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
#[8] "#666666"

## Create a new column named label to specify which genes should be labeled
df$label <- "NA"

## Set some criteria based on space in the figure and replace NA with gene names
df$label[df$diff.exp != "NO" & df$log2FoldChange < -2 & df$pval < 10^-50] <- rownames(df)[df$diff.exp != "NO" & df$log2FoldChange < -2 & df$pval < 10^-50]
df$label[df$diff.exp != "NO" & df$log2FoldChange > 2 & df$pval < 10^-50] <- rownames(df)[df$diff.exp != "NO" & df$log2FoldChange > 2 & df$pval < 10^-50]

## Plot
p <- ggplot(df, aes(x=log2FoldChange, y=-log10(pvalue), col=diff.exp)) + 
    geom_point() + xlim(-10,10) + #change xlim based on data
    scale_color_manual(values=c("#1B9E77","#666666","#D95F02")) +
    theme_test() +
    geom_text_repel(data=subset(df, label != "NA"), aes(label = label)) +#label only ones that are not NA
    labs(title = "GSK591 vs DMSO",
       x = "log2(fold change)",
       y = "-log10(P-value)")

pdf("volcano-plot.pdf", height=5, width=6)
p
dev.off()
