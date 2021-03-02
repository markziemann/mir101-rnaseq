up<-read.table("up.txt")
dn<-read.table("dn.txt")
boxplot(up$V1,dn$V1,log="y",names=c("up genes","dn genes"),ylab = "set annotations per gene")
mtext(paste("up mean=",signif(mean(up$V1),3),"dn mean=",signif(mean(dn$V1),3)))
wilcox.test(up$V1,dn$V1)

