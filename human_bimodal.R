x <- read.csv("GSE76381_EmbryoMoleculeCounts.cef.txt.gz",sep="\t",skip=1)
x <- x[-c(1:3),]
x <- x[,-2]
x[,-1] <- apply(x[,-1],2,as.numeric)
require(diptest)
P <- apply(apply(x[,-1],2,as.numeric),1,function(x){dip.test(x)$p.value})
table(p.adjust(P,"BH")<0.01)
D <- apply(apply(x[,-1],2,as.numeric),1,function(x){dip.test(x)$statistic})
x_human_bimodal <- x[order(D,decreasing=T)[1:200],]
#x_human_bimodal includes top ranked 200 bimodal genes for human. First column is the gene names.
