x <- read.csv("GSE76381_MouseEmbryoMoleculeCounts.cef.txt.gz",sep="\t",skip=2)
x <- x[-c(1:4),] #mouse
x <- x[,-2]
x[,-1] <- apply(x[,-1],2,as.numeric)
pca <- prcomp(scale(x[,-1]))
P <- pchisq(rowSums(scale(pca$x[,1:3])^2),2,lower.tail=F)
table(p.adjust(P,"BH")<0.01)
x_mouse <- x[p.adjust(P,"BH")<0.01,]
#x_mouse includes gene expression of selected genes. The first column is gene names
