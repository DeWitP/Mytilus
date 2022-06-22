source('/Users/Pierre/Desktop/Med19/May2021/EEMS_cluster2/plot_admixture_v5_function.R')
library(RcppCNPy)

# assembling the input table
dir="/Users/Pierre/Desktop/Med19/May2021/EEMS_cluster2"
pops="/inds2pops.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
#------------
setwd(dir)
tbl=data.frame(npyLoad("pcangsd.admix.Q.npy"))
npops=ncol(tbl)
i2p=read.table(paste(dir,pops,sep=""),header=F)

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)

ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance",vshift=0.05,angle=45,cex=0.8)

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
save(cluster.admix,file=paste("Med19-cluster2","_clusters.RData",sep=""))


# Trying to rearrange order of individuals, to look nicer:
#But of course the plotting function needs to be rewritten...

Qsorted <- read.delim(file="admix_data_table_Qsorted.txt", header=TRUE)
row.names(Qsorted)=Qsorted$ind
npops=2
vshift=0.05
hshift=0
angle=45
cex=0.8
require(colorspace)
require(RColorBrewer)
colors=c("tomato", "lightblue", "wheat","olivedrab", "cyan3","hotpink","gold","orange")
p=levels(Qsorted$pop)[1]
midpts=barplot(t(as.matrix(Qsorted[1:npops])), col=colors,xlab="", space=0,ylab="Ancestry", border=NA, xaxt="n",mgp=c(2,1,0))
pops=unique(as.character(Qsorted$pop))
np=0;lp=0;lpoints=c()
abline(v=0)
abline(h=1,lwd=1)
abline(h=0,lwd=1)
for( p in pops) {
  np0=table(Qsorted$pop==p)[2]
  lp0=np0/2
  lp=np+lp0
  np=np+np0
  lpoints=append(lpoints,lp)
  abline(v=np)
}
text(as.numeric(lpoints)+3.5+hshift,par("usr")[3]-0.1-vshift,labels=pops,xpd=TRUE, srt=angle, pos=2, cex=0.8)
row.names(Qsorted)=Qsorted$ind

