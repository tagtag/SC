x <- read.csv("GSE76381_EmbryoMoleculeCounts.cef.txt.gz",sep="\t",skip=1)
x <- x[-c(1:3),]
x <- x[,-2]
x[,-1] <- apply(x[,-1],2,as.numeric)
require(locfit)
Var <- apply(x[,-1],1,var)
Mean <- rowMeans(x[,-1])
Y <- Var^0.5/Mean
NLS <- nls(log10(Y)~0.5*log10(a/Mean+b),start=list(a=1,b=1),trace=T)
P <- pchisq(scale(log10(Y) - predict(NLS,Mean))^2,1,lower.tail=F) 
x_human_high <- x[!is.na(Y),][p.adjust(P,"BH")<0.01,]
#x_human_high includes highly variables genes for human. The first column is gene names.
