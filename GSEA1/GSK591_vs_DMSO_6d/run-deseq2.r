library('DESeq2')

sampleTable<-read.table("counts.txt", header=TRUE, check.names=F, row.names = 1)
head(sampleTable)

## DESeq2 takes whichever comes alphabetically first as the control group
samples <- data.frame(row.names=c("DMSO_rep1", "DMSO_rep2", "GSK_rep1", "GSK_rep2"), 
		      condition=as.factor(c(rep("DMSO",2),rep("GSK",2))))
samples

ddscounts <- DESeqDataSetFromMatrix(countData = sampleTable, colData=samples, design=~condition)
ddscounts
dds = DESeq(ddscounts)
res <- results(dds)
res<-res[order(res$padj),]
head(res)
f <-res[which(res$padj<0.05),]
write.table(file='diff-exp-fdr-0.05.txt', f, quote=F, col.names=T, sep='\t')
f <-res[which(res$padj<0.1),]
write.table(file='diff-exp-fdr-0.1.txt', f, quote=F, col.names=T, sep='\t')
f <-res[which(res$padj<0.01),]
write.table(file='diff-exp-fdr-0.01.txt', f, quote=F, col.names=T, sep='\t')
f <-res[which(res$padj<=1),]
write.table(file='diff-exp-fdr.txt', f, quote=F, col.names=T, sep='\t')

png(filename="deseq2_MAplot.png")
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.off()

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

pdf('deseq2_dispersion.pdf', height=4, width=4)
plotDispEsts(dds,main="Dispersion")
dev.off()

pdf('deseq2_pca.pdf', height=6, width=9)
plotPCA(rld, intgroup=c("condition"))
dev.off()

sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
library("gplots")
library("RColorBrewer")
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf('deseq2_sampledistance.pdf', pointsize=8)
heatmap.2(sampleDistMatrix, trace="none", col=colours, margins =c(12,12))
dev.off()

norm.counts <-as.data.frame(counts(dds, normalized=T))
write.table(norm.counts, file='norm-counts-all-data.txt', quote=F, col.names=T, sep='\t')
