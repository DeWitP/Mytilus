# installing required packages: do this once, then remark with #
#install.packages("vegan")
#install.packages("ggplot2")
# install.packages("RcppCNPy")
# install.packages("pheatmap")
library(RcppCNPy)

setwd("/Users/Pierre/Desktop/Med19/May2021/EEMS_cluster2") # change this to where your scp'd files are
bams=read.table("bams.nc")[,1] # list of bam files
bams0=bams
# removing leading / trailing filename bits from sample names
bams=sub(".fastq.trim.bam","",bams)
bams=sub(".+/","",bams)
length(bams)

#------ PCAngsd magic

# kinship (from pcangsd)
kin= npyLoad("pcangsd.kinship.npy") 
dimnames(kin)=list(bams, bams)
pheatmap::pheatmap(kin)
plot(hclust(as.dist(1-kin),method="ave"))

# inbreeding (from pcangsd)
inbr= npyLoad("pcangsd.inbreed.npy") 

# reading pcangsd covariance matrix, converting to correlation-based distances
pcc = as.matrix(read.table("pcangsd.cov"))
dimnames(pcc)=list(bams,bams)
pccd=1-cov2cor(pcc)

# must run admixturePlotting_v5.R first for the following to work
load('Med19-cluster2_clusters.RData')

# reading the ibs matrix
ma = as.matrix(read.table("OK.ibsMat"))
dimnames(ma)=list(bams,bams)

# population designations (first letter of sample names)
i2p=read.table(file="inds2pops.txt")
pops<-i2p[,2]

# ---- sanity checks

# heatmaps of distances
library(pheatmap)
pheatmap(ma)
pheatmap(kin)
pheatmap(pccd)
# which one looks more structured?

# hierarchical clustering trees
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.8) # clustering of samples by IBS 
hc=hclust(as.dist(pccd),"ave")
plot(hc,cex=0.8) # clustering of samples by SNP correlation (from pcangsd) 

#----- vegan stuff (on IBS (ma) or pcangsd-derived (pcccd) distances)

library(vegan)

# PCoA
#ordi=capscale(ma~1)
ordi=capscale(pccd~1)

# eigenvalue
plot(ordi$CA$eig) 
# how many "interesting" ordination axes do we see?

# plotting vegan way (I really like ordispider and ordiellipse functions)

# extracting "scores" table, to plot
axes2plot=c(1) # which MDS to plot
scores=data.frame(scores(ordi,display="sites",choices=axes2plot))

# making color scheme for our groups (colramp trick makes it work with >12 groups)
library(RColorBrewer)
#groups2color=cluster.admix
groups2color=pops
colramp = colorRampPalette(brewer.pal(length(unique(groups2color)),"Set2"))
myColors=colramp(length(unique(groups2color)))
names(myColors) = sort(unique(groups2color))

plot(scores$MDS1~groups2color,pch=16,col="grey",asp=1, ylab="MDS1", ylim=c(-3,3), xaxt="n")
## Draw x-axis without labels.
axis(side = 1, labels = FALSE)

## Draw the x-axis labels.
text(x = 1:length(myColors),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = names(myColors),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Decrease label size.
     cex = 0.8)

# importing and aligning coverage data (from quality.txt file produced by plotQC.R)
cover=read.table("quality.txt")
cover$V1=sub(".fastq.trim.bam","",cover$V1)
row.names(cover)=sub(".+/","",cover$V1)
cover$V1=NULL
cover=cover[bams,]

#--------- statistical tests (PERMANOVA)

# does coverage and population designation aligns with our "interesting" axes? (use this function to check is any parameter aligns with ordination) 

# assembling data frame of variables to corelate with ordination (cover and pops)
pc=data.frame(cbind(pops,cover),stringsAsFactors=F)
pc$cover=as.numeric(pc$cover)

# fitting them to chosen ordination axes (let's take the first 3)
ef=envfit(ordi,pc,choices=c(1:3))
ef
# note that pop loads on MDS1, cover - on MDS3 (can you make a plot to verify that?). 


# How much *overall* variation is attributable to sampling site and coverage? Are they significant factors?
adonis(pccd~pops,pc)
# factor "pops" explains 0.574% of variation and is significant. 
# factor "cover" explains 0.377% of variation and is significant..
 
# -------- can we remove effect of coverage from ordination?

# partial ordination - removing effect of coverage
# pp1=capscale(ma~1+Condition(cover))
pp1=capscale(pccd~1+Condition(cover))

# let's plot this like before
axes2plot=c(1) 
scores=data.frame(scores(pp1,display="sites",choices=axes2plot))

plot(scores$MDS1~groups2color,pch=16,col="grey",asp=1, ylab="MDS1", ylim=c(-3,3), xaxt="n")
## Draw x-axis without labels.
axis(side = 1, labels = FALSE)

## Draw the x-axis labels.
text(x = 1:length(myColors),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = names(myColors),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Decrease label size.
     cex = 0.8)

# Can we do a Tukey posthoc?
library(car)
model<-glm(scores$MDS1~groups2color)
test<-aov(model)
TUK <- TukeyHSD(test, ordered=TRUE, conf.level=0.95)
