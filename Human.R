x <- read.csv("GSE76381_EmbryoMoleculeCounts.cef.txt.gz",sep="\t",skip=1)
x <- x[-c(1:3),]
x <- x[,-2]
x[,-1] <- apply(x[,-1],2,as.numeric)
pca <- prcomp(scale(x[,-1]))
P <- pchisq(rowSums(scale(pca$x[,1:2])^2),2,lower.tail=F)
table(p.adjust(P,"BH")<0.01)
x_human <- x[p.adjust(P,"BH")<0.01,]
#x_human includes expression profiles of selected genes. The fiest column is gene names
