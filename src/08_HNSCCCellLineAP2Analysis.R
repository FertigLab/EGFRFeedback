library('ProjectTemplate')
load.project()


# HNSCC gene expression
HNSCC.Gene <- t(apply(HNSCC.CellLines.FRMA[geneProbeSelect,
                                           HNSCC.CellLines.sampleAnnot$CelFile],1,function(x)
                                           {tapply(x,HNSCC.CellLines.sampleAnnot$CellLine_Sen,mean)}))
row.names(HNSCC.Gene) <- names(geneProbeSelect)
colnames(HNSCC.Gene) <- sapply(strsplit(colnames(HNSCC.Gene),split="_"),
                               function(x){paste0(x[1:(length(x)-1)],collapse="_")})


# cell lines with pre and post-treatment information
CTX.Names <- grep('CTX',colnames(HNSCC.Gene),value=T)
CTX.Names <- CTX.Names[grep('1CC8',CTX.Names,invert=T)]

CTX <- tapply(HNSCC.CellLines.sampleAnnot$Hiro.CTX,HNSCC.CellLines.sampleAnnot$CellLine,unique)


# HPV-positive and HPV-negative cells
HPVNeg <- c('SCC1','SCC25','SQ20B','SCC61')
HPVNeg <- names(sort(CTX[HPVNeg]))

HPVPos <- c('SCC047', 'SCC090', '93VU147T')
HPVPos <- names(sort(CTX[HPVPos]))


# compute the average EGFR family expression

EGFRExprs <- matrix(NA, nrow=length(CTX.Names), ncol=2,
                    dimnames=list(sapply(strsplit(CTX.Names,split="_"),function(x){x[[1]]}),
                                  c('Control', 'Cetuximab')))
EGFRMin <- EGFRMax <- EGFRExprs

for (c in row.names(EGFRExprs)) {
  for (e in colnames(EGFRExprs)) {
    samp <- colnames(HNSCC.CellLines.FRMA)
    samp <- samp[HNSCC.CellLines.sampleAnnot[samp,'CellLine'] == c]
    samp <- samp[HNSCC.CellLines.sampleAnnot[samp,'Treatment..0.No.Tx..1.Cetux.Tx.24.hr.']==
                    as.numeric(e == 'Cetuximab')]
    EGFRExprs[c,e] <- mean(HNSCC.CellLines.FRMA[geneProbeSelect['EGFR'],samp])
    EGFRMin[c,e] <- min(HNSCC.CellLines.FRMA[geneProbeSelect['EGFR'],samp])
    EGFRMax[c,e] <- max(HNSCC.CellLines.FRMA[geneProbeSelect['EGFR'],samp])    
  }
}

delEGFR <- sweep(EGFRExprs,1,EGFRExprs[,'Control'])


# changes in EGFR in HNSCC cell lines
pdf('graphs/DelEGFRHNSCC.pdf')
plot(CTX[row.names(delEGFR)], 
     delEGFR[,'Cetuximab'], 
     col=ifelse(row.names(delEGFR) %in% HPVPos, 'blue','red'), 
     pch=ifelse(row.names(delEGFR) %in% HPVPos, 15,19),
     ylim=c(-1,1),xlim=c(0,130),
     ylab='EGFR exprs in Trt - EGFR exprs in Control',
     xlab='Survival (%)')
text(CTX[row.names(delEGFR)], 
     delEGFR[,'Cetuximab'],row.names(delEGFR),pos=2,
     col=ifelse(row.names(delEGFR) %in% HPVPos, 'blue','red'))

abline(lm(as.numeric(delEGFR[HPVNeg,'Cetuximab'])
          ~ as.numeric(CTX[HPVNeg])),
       lty=2, lwd=2, col='red')

abline(lm(as.numeric(delEGFR[HPVPos,'Cetuximab'])
          ~ as.numeric(CTX[HPVPos])),
       lty=2, lwd=2, col='blue')

title(sprintf('HPV- %0.2f, HPV+ %0.2f',
              cor.test(as.numeric(CTX[HPVNeg]),as.numeric(delEGFR[HPVNeg,'Cetuximab']),
                       alternative='less')$p.value,
              cor.test(as.numeric(CTX[HPVPos]),as.numeric(delEGFR[HPVPos,'Cetuximab']),
                       alternative='less')$p.value))

legend('topright', col=c('blue','red'),pch=c(15,19),
       c('HPV+','HPV-'))
dev.off()




# Gene set statistics for AP-2 for HNSCC cell lines
load('data/TRANSFAC_Genes_2014.Rda')

# analysis across panel of cell lines
AP2Up <- read.table('reports/LINCS2000/AP2Up.grp',header=T)
AP2Down <- read.table('reports/LINCS2000/AP2Down.grp',header=T)

#delEGFRHNSCC <- HNSCC.Gene['EGFR',paste0(c(HPVNeg,HPVPos),'_CTX24h')] - 
#  HNSCC.Gene['EGFR',c(HPVNeg,HPVPos)]

#delAllHNSCC <- HNSCC.Gene[,paste0(c(HPVNeg,HPVPos),'_CTX24h')] - 
#  HNSCC.Gene[,c(HPVNeg,HPVPos)]

## HPV-negative
delEGFRHNSCC <- HNSCC.Gene['EGFR',paste0(c(HPVNeg),'_CTX24h')] - 
  HNSCC.Gene['EGFR',c(HPVNeg)]

delAllHNSCC <- HNSCC.Gene[,paste0(c(HPVNeg),'_CTX24h')] - 
  HNSCC.Gene[,c(HPVNeg)]

delAllHNSCC <- delAllHNSCC[setdiff(row.names(delAllHNSCC),'EGFR'),]

mm <- model.matrix(~delEGFRHNSCC)
delFit <- eBayes(lmFit(delAllHNSCC[,names(delEGFRHNSCC)],mm))

AP2StatHNSCCUp <- geneSetTest(statistics = delFit$t[,2],
                            index=intersect(AP2Up[,1],
                                            row.names(delAllHNSCC)),
                            alternative='greater')
AP2StatHNSCCDown <- geneSetTest(statistics = delFit$t[,2],
                                index=intersect(AP2Down[,1],
                                                row.names(delAllHNSCC)),
                                alternative='less')


pdf('graphs/HNSCCCellLinesEnrichment.pdf')
par(mfrow=c(2,1))
plot(sort(delFit$t[,2]),type='l',ylab='t EGFR association',ylim=c(-10,10))
  
for (p in setdiff(AP2Up[,1],'EGFR')) {
  points(x=which(names(sort(delFit$t[,2]))  ==p),
         y=delFit$t[p,2],
         pch=19,col='red')
  text(x=which(names(sort(delFit$t[,2]))  ==p),
       y=delFit$t[p,2]+0.5,
       labels=p,col='red',srt = 90,pos=4,offset=0)
}
  
for (p in setdiff(AP2Down[,1],'EGFR')) {
  points(x=which(names(sort(delFit$t[,2]))  ==p),
         y=delFit$t[p,2],
         pch=19,col='blue')
  text(x=which(names(sort(delFit$t[,2]))  ==p),
       y=delFit$t[p,2]-0.5,
       labels=p,col='blue',srt = 90,pos=2,offset=0)
}

title(sprintf('HPV-negative Up %0.2f, Down %0.2f',
              AP2StatHNSCCUp,AP2StatHNSCCDown))
  
## HPV-positive

delEGFRHNSCC <- HNSCC.Gene['EGFR',paste0(c(HPVPos),'_CTX24h')] - 
  HNSCC.Gene['EGFR',c(HPVPos)]

delAllHNSCC <- HNSCC.Gene[,paste0(c(HPVPos),'_CTX24h')] - 
  HNSCC.Gene[,c(HPVPos)]

delAllHNSCC <- delAllHNSCC[setdiff(row.names(delAllHNSCC),'EGFR'),]

mm <- model.matrix(~delEGFRHNSCC)
delFit <- eBayes(lmFit(delAllHNSCC[,names(delEGFRHNSCC)],mm))

AP2StatHNSCCUp <- geneSetTest(statistics = delFit$t[,2],
                              index=intersect(AP2Up[,1],
                                              row.names(delAllHNSCC)),
                              alternative='greater')
AP2StatHNSCCDown <- geneSetTest(statistics = delFit$t[,2],
                                index=intersect(AP2Down[,1],
                                                row.names(delAllHNSCC)),
                                alternative='less')



plot(sort(delFit$t[,2]),type='l',ylab='t EGFR association',ylim=c(-10,10))

for (p in setdiff(AP2Up[,1],'EGFR')) {
  points(x=which(names(sort(delFit$t[,2]))  ==p),
         y=delFit$t[p,2],
         pch=19,col='red')
  text(x=which(names(sort(delFit$t[,2]))  ==p),
       y=delFit$t[p,2]+0.5,
       labels=p,col='red',srt = 90,pos=4,offset=0)
}

for (p in setdiff(AP2Down[,1],'EGFR')) {
  points(x=which(names(sort(delFit$t[,2]))  ==p),
         y=delFit$t[p,2],
         pch=19,col='blue')
  text(x=which(names(sort(delFit$t[,2]))  ==p),
       y=delFit$t[p,2]-0.5,
       labels=p,col='blue',srt = 90,pos=2,offset=0)
}

title(sprintf('HPV-positive Up %0.2f, Down %0.2f',
              AP2StatHNSCCUp,AP2StatHNSCCDown))
dev.off()


### heatmap
delEGFRHNSCC <- HNSCC.Gene['EGFR',paste0(c(HPVNeg,HPVPos),'_CTX24h')] - 
  HNSCC.Gene['EGFR',c(HPVNeg,HPVPos)]

delAllHNSCC <- HNSCC.Gene[,paste0(c(HPVNeg,HPVPos),'_CTX24h')] - 
  HNSCC.Gene[,c(HPVNeg,HPVPos)]

delAllHNSCC <- delAllHNSCC[setdiff(row.names(delAllHNSCC),'EGFR'),]

mm <- model.matrix(~delEGFRHNSCC)
delFit <- eBayes(lmFit(delAllHNSCC[,names(delEGFRHNSCC)],mm))

AP2StatHNSCCUp <- geneSetTest(statistics = delFit$t[,2],
                              index=intersect(AP2Up[,1],
                                              row.names(delAllHNSCC)),
                              alternative='greater')
AP2StatHNSCCDown <- geneSetTest(statistics = delFit$t[,2],
                                index=intersect(AP2Down[,1],
                                                row.names(delAllHNSCC)),
                                alternative='less')

delPlot <- rbind(EGFR=delEGFRHNSCC,
                 delAllHNSCC[intersect(TF2Gene[['AP-2alpha']],
                                       row.names(delAllHNSCC)),])
pdf('graphs/DelEGFRHNSCCHeatmap.pdf')

heatmap.2(delPlot[,paste0(c(HPVNeg,HPVPos),'_CTX24h')],
          scale='row',trace='none',col=greenred, Colv=F,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          labCol = c(HPVNeg,HPVPos),
          ColSideColors = rep(c('red','blue'),
                              c(length(HPVNeg),
                                length(HPVPos))),
          colRow=ifelse(row.names(delPlot) %in% AP2Up[,1],'red','black'),
          main=sprintf('Up p %0.2f, Down p %0.2f',
                       AP2StatHNSCCUp, AP2StatHNSCCDown))
dev.off()

