library('ProjectTemplate')
load.project()

# compare between replicate variablity to between condition variablity
# to select probes for analysis
Premalignancy.control <- 
  row.names(Premalignancy.sampleAnnot)[Premalignancy.sampleAnnot$HAND.ID == 'FaDu']
HNSCC.CellLines.control <- row.names(HNSCC.CellLines.sampleAnnot)[
  HNSCC.CellLines.sampleAnnot$CellLine %in% 
    c(unique(HaCaT.sampleAnnot$CellType),'FaDu') &
    HNSCC.CellLines.sampleAnnot$Treatment..0.No.Tx..1.Cetux.Tx.24.hr.==0]
HaCaT.control <- row.names(HaCaT.sampleAnnot)[HaCaT.sampleAnnot$CellLine %in%
                                                HNSCC.CellLines.sampleAnnot$CellLine]

control.FRMA <- cbind(Premalignancy.FRMA[,Premalignancy.control],
                      HNSCC.CellLines.FRMA[,HNSCC.CellLines.control],
                      HaCaT.FRMA[,HaCaT.control])
colnames(control.FRMA) <- c(Premalignancy.control,
                            HNSCC.CellLines.control,
                            HaCaT.control)

control.sampleAnnot <- c(as.character(Premalignancy.sampleAnnot[
  Premalignancy.control,'HAND.ID']),
                         HNSCC.CellLines.sampleAnnot[HNSCC.CellLines.control,
                                                     'CellLine'],
                         HaCaT.sampleAnnot[HaCaT.control,'CellLine'])
names(control.sampleAnnot) <- c(Premalignancy.control,
                                HNSCC.CellLines.control,
                                HaCaT.control)

control.FRMA.diff <- control.FRMA
control.mean <- matrix(NA, nrow=nrow(control.FRMA), 
                       ncol=length(unique(control.sampleAnnot)),
                       dimnames=list(row.names(control.FRMA),
                                     unique(control.sampleAnnot)))
for (s in colnames(control.mean)) {
  sIdx <- names(control.sampleAnnot)[control.sampleAnnot %in% s]
  control.mean[,s] <- apply(control.FRMA[,sIdx],1,mean)
  control.FRMA.diff[,sIdx] <- sweep(control.FRMA.diff[,sIdx],1,
                                    control.mean[,s])
}

HaCaT.mean <- matrix(NA, nrow=nrow(control.FRMA),
                     ncol=length(unique(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),'ExperimentalConditions'])),
                     dimnames=list(row.names(control.FRMA),
                                   unique(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),'ExperimentalConditions'])))
HaCaT.diff <- HaCaT.FRMA.pSVA
for (s in colnames(HaCaT.mean)) {
  sIdx <- colnames(HaCaT.FRMA.pSVA)[HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),'ExperimentalConditions'] == s]
  if (length(sIdx) > 1) {
    HaCaT.mean[,s] <- apply(HaCaT.FRMA.pSVA[,sIdx],1,mean)
    HaCaT.diff[,sIdx] <- sweep(HaCaT.FRMA.pSVA[,sIdx],1, HaCaT.mean[,s])
  } else {
    HaCaT.mean[,s] <- HaCaT.FRMA.pSVA[,sIdx]
    HaCaT.diff[,sIdx] <- NA
  }
}

exp.mad <- apply(cbind(HaCaT.mean,control.mean),1,mad)
rep.mad <- apply(cbind(HaCaT.diff,control.FRMA.diff),1,mad,na.rm=T)

pdf('graphs/Supplement/MADProbesRep.pdf')
plot(exp.mad,rep.mad, 
     xlab='MAD of experimental conditions (log2 exprs)',
     ylab='MAD of replicates (log2 exprs)')
dev.off()

DevFN <- function(x){return(exp.mad[x]/(rep.mad[x] + .Machine$double.eps))}

genes2Probes <- revmap(as.list(hgu133plus2SYMBOL))

geneProbeSelect <- sapply(genes2Probes,
                          function(x){names(which.max(DevFN(x)))})

HaCaT.Gene.FRMA.pSVA <- HaCaT.FRMA.pSVA[geneProbeSelect,]
row.names(HaCaT.Gene.FRMA.pSVA) <- names(geneProbeSelect)

HaCaT.Gene.FRMA.pSVA <- HaCaT.Gene.FRMA.pSVA[substr(row.names(HaCaT.Gene.FRMA.pSVA),1,3)!='LOC',]
HaCaT.Gene.FRMA.pSVA <- HaCaT.Gene.FRMA.pSVA[grep('orf',row.names(HaCaT.Gene.FRMA.pSVA),value=T,invert=T),]


pdf('graphs/clusterGenes.pdf')
for (c in colnames(HaCaT.sampleAnnot)) {
  plotColoredClusters(standard.pearson(HaCaT.Gene.FRMA.pSVA),
                      labs=HaCaT.sampleAnnot[colnames(HaCaT.Gene.FRMA.pSVA),'ExperimentalConditions'],
                      cols=as.character(as.numeric(as.factor(HaCaT.sampleAnnot[colnames(HaCaT.Gene.FRMA.pSVA),c]))), 
                      xlab='')
  if (all(is.na(HaCaT.sampleAnnot[colnames(HaCaT.Gene.FRMA.pSVA),c]))) next
  legend('topleft', fill=as.character(1:length(unique(HaCaT.sampleAnnot[colnames(HaCaT.Gene.FRMA.pSVA),c]))), 
         title=c,
         legend=levels(factor(HaCaT.sampleAnnot[colnames(HaCaT.Gene.FRMA.pSVA),c])))
}
dev.off()

# cache the gene level data
ProjectTemplate::cache('HaCaT.Gene.FRMA.pSVA')
ProjectTemplate::cache('geneProbeSelect')

# output the gene probe selection as a supplemental table
geneProbeSelectOut <- cbind(probe=geneProbeSelect,gene=names(geneProbeSelect))
write.table(geneProbeSelectOut,file='reports/GeneProbeSelection.txt',sep="\t",row.names=F,quote=F)
