path = "."

infiles <- dir(pattern='\\.xls$')


plot.GSEA <- function(file){
  x <-read.table(file,header=T,row.names=1)

  y<-head(x[order(x$ES),],n=20L)
  y<-y[which(y$ES<0),]

#newdata <- mydata[ which(gender=='F' & age > 65),]
#detach(newdata)

  z<-head(x[order(-x$ES),],n=20L)
  z<-z[which(z$ES>0),]

  df <- rbind(y,z)
  df<-df[order(df$ES),]
  barplot(df$ES,main=file,horiz=TRUE,names.arg=row.names(df))
}

pdf(file="gsea_results.pdf",width=15,height=10)
par(las=2) ; par(mar=c(10,50,1,1))
lapply(infiles,plot.GSEA)
dev.off()

q()

path = "."
infiles <- dir(pattern='\\.xls$')

plot.GSEA <- function(file){
  x <-read.table(file,header=T,row.names=1)
  y<-head(x[order(x$ES),],n=20L)
  z<-head(x[order(-x$ES),],n=20L)
  df <- rbind(y,z)
  df<-df[order(df$ES),]
  pdf(file=sub("\\.xls$",".pdf", file),width=15,height=10)
  par(las=2) ; par(mar=c(10,50,1,1))
  barplot(df$ES,main=c(file=sub("\\.xls$",""),horiz=TRUE,names.arg=rev(row.names(df)))
}

lapply(infiles,plot.GSEA)

q()
#
x<-read.table("pairedeffectof149inHG.xls.fmt.rnk.c2.cp.reactome.v5.1.symbols.gmt.GseaPreranked.1457314982486.xls",header=T,row.names=1)
y<-head(x[order(x$ES),],n=20L)
z<-head(x[order(-x$ES),],n=20L)
df <- rbind(y,z) 
df<-df[order(df$ES),]
pdf(file="barplot.pdf",width=15,height=10) ; par(las=2) ; par(mar=c(10,50,1,1)) ; barplot(df$ES,main="Enrichment Score",horiz=TRUE,names.arg=rev(row.names(df))) ; dev.off()

