x <- read.csv("GSE76381_EmbryoMoleculeCounts.cef.txt.gz",sep="\t",skip=1)
x <- x[-c(1:3),]
x <- x[,-2]
x[,-1] <- apply(x[,-1],2,as.numeric)
X <- x[,-1]
colnames(X)<-colnames(x[,-1])
rownames(X)<-x[,1]
require(monocle)
pd <- new("AnnotatedDataFrame", data = data.frame(t(X)))
fd <- new("AnnotatedDataFrame", data = data.frame(X))
HSMM <- newCellDataSet(as.matrix(X),phenoData = pd, featureData = fd)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
fData(HSMM)$use_for_ordering <-fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- reduceDimension(HSMM,
                              max_components = 2,
                              norm_method = 'log',
                              num_dim = 3,
                              reduction_method = 'tSNE',
                              verbose = T)
HSMM <- clusterCells(HSMM,
                 rho_threshold = 2,
                 delta_threshold = 4,
                 skip_rho_sigma = T,
                 verbose = F)
HSMM_expressed_genes <-  row.names(subset(fData(HSMM),num_cells_expressed >= 10))
clustering_DEG_genes <-
    differentialGeneTest(HSMM[HSMM_expressed_genes,],
          fullModelFormulaStr = '~Cluster',
          cores = 1)
table(clustering_DEG_genes$qval<0.01)
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]
#HSMM_ordering_genes includes top ranked 200 genes that dpFature selected.
