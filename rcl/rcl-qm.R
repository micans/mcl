
library(Matrix)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(dplyr))

# quickMarkers was copied verbatim from Matthew Young's SoupX package:
# https://github.com/constantAmateur/SoupX/

quickMarkers = function(toc,clusters,N=10,FDR=0.01,expressCut=0.9){
  #Convert to the more manipulable format
  toc = as(toc,'dgTMatrix')
  w = which(toc@x>expressCut)
  #Get the counts in each cluster
  clCnts = table(clusters)
  nObs = split(factor(rownames(toc))[toc@i[w]+1],clusters[toc@j[w]+1])
  nObs = sapply(nObs,table)
  #Calculate the observed and total frequency
  nTot = rowSums(nObs)
  tf = t(t(nObs)/as.integer(clCnts[colnames(nObs)]))
  ntf = t(t(nTot - nObs)/as.integer(ncol(toc)-clCnts[colnames(nObs)]))
  idf = log(ncol(toc)/nTot)
  score = tf*idf
  #Calculate p-values
  qvals = lapply(seq_len(ncol(nObs)),function(e)
                 p.adjust(phyper(nObs[,e]-1,nTot,ncol(toc)-nTot,clCnts[colnames(nObs)[e]],lower.tail=FALSE),method='BH'))
  qvals = do.call(cbind,qvals)
  colnames(qvals) = colnames(nObs)
  #Get gene frequency of second best
  sndBest = lapply(seq_len(ncol(tf)),function(e) apply(tf[,-e,drop=FALSE],1,max))
  sndBest = do.call(cbind,sndBest)
  colnames(sndBest) = colnames(tf)
  #And the name
  sndBestName = lapply(seq_len(ncol(tf)),function(e) (colnames(tf)[-e])[apply(tf[,-e,drop=FALSE],1,which.max)])
  sndBestName = do.call(cbind,sndBestName)
  colnames(sndBestName) = colnames(tf)
  rownames(sndBestName) = rownames(tf)
  #Now get the top N for each group
  w = lapply(seq_len(ncol(nObs)),function(e){
             o = order(score[,e],decreasing=TRUE)
             if(sum(qvals[,e]<FDR)>=N){
               o[seq(N)]
             }else{
               o[qvals[o,e]<FDR]
             }
                 })
  #Now construct the data.frame with everything
  ww = cbind(unlist(w,use.names=FALSE),rep(seq_len(ncol(nObs)),lengths(w)))
  out = data.frame(gene = rownames(nObs)[ww[,1]],
                   cluster = colnames(nObs)[ww[,2]],
                   geneFrequency = tf[ww],
                   geneFrequencyOutsideCluster = ntf[ww],
                   geneFrequencySecondBest = sndBest[ww],
                   geneFrequencyGlobal = nTot[ww[,1]]/ncol(toc),
                   secondBestClusterName = sndBestName[ww],
                   tfidf = score[ww],
                   idf = idf[ww[,1]],
                   qval = qvals[ww],
                   stringsAsFactors=FALSE)
  return(out)
}


args <- R.utils::commandArgs(trailingOnly=TRUE)

  ##
  ## This script is created for use by rcl-qc.sh, so unless further incentives
  ## arrive, eight positional parameters remain the amazing interface.
  ## 

if (length(args) != 8) {
  stop("Need <MATRIXFILE> <GENEROWNAMES> <CLUSTERFILE> <OUTPUTFILE> <QMCACHEFILE> <QMCACHEDATA> <NUM> <FDR>")
}

matrixfile   <- args[[1]]
generownames <- args[[2]]
clusterfile  <- args[[3]]
outputfile   <- args[[4]]
qmcacheqm    <- args[[5]]         # if present read, otherwise create.
qmcachedata  <- args[[6]]         # if present read, otherwise create.
qmN          <- as.numeric(args[[7]])
qmFDR        <- as.numeric(args[[8]])

cls <- read.table(clusterfile, row.names=1, sep="\t")
geneNames <- readLines(generownames)


if (!file.exists(qmcachedata) || !file.exists(qmcacheqm)) {

  cat("Reading matrix ..", file=stderr())
    themtx <- readMM(matrixfile)
    dimnames(themtx) <- list(geneNames, rownames(cls))
  cat(".. done\n", file=stderr())

  cat("Running permissive quickMarkers with N=5 FDR=0.01 ..", file=stderr())
    qm <- quickMarkers(themtx, cls$V2, N=5, FDR=0.01)
  cat(".. done\n", file=stderr())
  cat(sprintf("Writing qm and data cache files %s, %s\n", qmcacheqm, qmcachedata), file=stderr())
    write.table(qm, file=qmcacheqm, sep="\t", quote=FALSE)
    u  <- unique(qm$gene)
    themtx <- t(as.matrix(themtx[u,]))
    write.table(as.matrix(themtx), file=qmcachedata, sep="\t", quote=FALSE)

} else {
  cat(sprintf("Reading qm and data cache files %s, %s ..", qmcacheqm, qmcachedata), file=stderr())
    # _ Note: not a sparse matrix.
    themtx <- read.table(qmcachedata, header=T, row.names = 1, check.names=FALSE, as.is=TRUE)
    qm <- read.table(qmcacheqm, header=T, check.names=FALSE, as.is=TRUE)
  cat(".. done\n", file=stderr())
}

# cat(sprintf("Selecting quickMarkers with N=%d FDR=%f\n", qmN, qmFDR), file=stderr())
# qm2 <- as.data.frame(qm %>% group_by(cluster) %>% slice_min(order_by = qval, n = qmN))
# qm2 <- qm2[qm2$qval <= qmFDR,]

cat(sprintf("Selecting quickMarkers with N=%d tfidf=%f\n", qmN, qmFDR), file=stderr())
qm2 <- as.data.frame(qm %>% group_by(cluster) %>% slice_max(order_by = tfidf, n = qmN))
qm2 <- qm2[qm2$tfidf >= qmFDR,]

write.table(as.matrix(themtx[,unique(qm2$gene)]), file=outputfile, sep="\t", quote=FALSE)


