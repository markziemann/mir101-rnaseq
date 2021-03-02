#run deseq2 with Matt's counts from AGRF
library("tidyverse")
library("DESeq2")
library("fgsea")
library("gplots")

#read counts
counts<-read.table("CAGRF13809_cnt.tsv",header=T,row.names=1)

#curate sample sheet. miR101=MM1 Control=MMN
des<-as.data.frame(colnames(counts))
des$grp<-as.numeric(grepl("MM1",colnames(counts)))
rownames(des)<-des[,1]
des[,1]=NULL

dds <- DESeqDataSetFromMatrix(countData = counts , colData = des, design = ~ grp)
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)

#stick on the normalised expression values to the table
zz<-cbind(as.data.frame(z),assay(vsd))

#sort by p-value
mm1<-as.data.frame(zz[order(zz$pvalue),])

#some plots
pdf("CAGRF13809_plots.pdf")
sig<-subset(zz,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("Ctrl vs miR-101:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")
plot(log2(zz$baseMean),zz$log2FoldChange,cex=0.6, xlab="log2 base mean", ylab="log2 fold change" ,pch=19,col="#838383")
points(log2(sig$baseMean),sig$log2FoldChange,cex=0.6,pch=19,col="red")
mtext(HEADER)

top<-head(sig,20)
text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)

#volcano plot
plot(zz$log2FoldChange, -log2(zz$pvalue) ,cex=0.6, xlim=c(-4,6),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")
points(sig$log2FoldChange, -log2(sig$pvalue),cex=0.6,pch=19,col="red")
text(top$log2FoldChange+0.5, -log2(top$pvalue), labels = rownames(top),cex=0.7)
mtext(HEADER)

# top N genes
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(  as.matrix(zz[1:100,c(7:ncol(zz))]), col=colfunc(25),scale="row", 
 trace="none",margins = c(6,6), cexRow=.4, main="Top 100 genes")
dev.off()

#output DGE table
write.(zz,file="CAGRF13809_deseq.tsv",quote=F,sep="\t")

# score DE results by pvalue
rnk<-as.data.frame(mm1$log2FoldChange*-log(mm1$pvalue))
rownames(rnk)<-row.names(mm1)
colnames(rnk)="Score"
write.table(rnk,file="CAGRF13809.rnk",sep='\t',quote=F)

#now follow Stephen Turner example
ranks <- deframe(rnk)
names(ranks)<-rownames(rnk)
ranks2<-rank(ranks)-length(which(ranks<0))

#load in genesets
pathways.msigdb <- gmtPathways("../genesets/msigdb.v6.2.symbols.gmt")

# run fgsea
fgseaRes<-fgsea(pathways=pathways.msigdb,stats=ranks2,minSize=10,nperm=10000)
fgseaRes<-fgseaRes[order(fgseaRes$pval),]
fgseaRes<-subset(fgseaRes,padj<0.05)
fgseaRes<-head(fgseaRes,100)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
