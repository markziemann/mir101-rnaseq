#install.packages("devtools")
#library("devtools")
#devtools::install_github("markziemann/Mitch")
library("mitch")

# import rank files for data
mir101ox<-read.table("CAGRF13809.rnk",header=T,row.names=1)
colnames(mir101ox)<-"mir101ox"
mst1tg<-read.table("WTvsMSTTG.tsv.rnk.hum",header=T,row.names=1)
colnames(mst1tg)<-"mst1tg"
x<-merge(mir101ox,mst1tg,by=0)
rownames(x)<-x$Row.names
x$Row.names=NULL

# use reactomes
genesets<-gmt_import("ReactomePathways.gmt")

# run mitch with both confidence and significance
res<-mitch_calc(x,genesets,resrows=100,bootstraps=1000,priority="confidence",cores=6)
mitch_plots(res,outfile="mir101oxVSmst1tg_conf.pdf")
mitch_report(res,"mir101oxVSmst1tg_conf.html")
res_conf<-res

res<-mitch_calc(x,genesets,resrows=100,bootstraps=1000,priority="significance",cores=6)
mitch_plots(res,outfile="mir101oxVSmst1tg_sig.pdf")
mitch_report(res,"mir101oxVSmst1tg_sig.html")
res_sig<-res

save.image("mir101osVSmst1tg_mitch.RData")
