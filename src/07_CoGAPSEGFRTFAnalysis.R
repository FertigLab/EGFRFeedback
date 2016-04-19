library('ProjectTemplate')
load.project()

# get transcription factor targets for 
# TF analysis of each of the CoGAPS paterns
load('data/TRANSFAC_Genes_2014.Rda')

TF2Gene <- lapply(TF2Gene,intersect,row.names(CoGAPS_SDFIVE_3$Amean))
TF2Gene <- TF2Gene[sapply(TF2Gene,length)>=5]
TF2Gene <- TF2Gene[sapply(TF2Gene,length)<100]

# format TF2 gene list for output
TF2GeneTable <- matrix(F,nrow=nrow(CoGAPS_SDFIVE_3$Amean),
                       ncol=length(TF2Gene),
                       dimnames = list(row.names(CoGAPS_SDFIVE_3$Amean),
                                       names(TF2Gene)))
for (tf in names(TF2Gene)) {
  TF2GeneTable[TF2Gene[[tf]],tf] <- T
}

write.table(data.frame(Gene=row.names(TF2GeneTable),
                       TF2GeneTable),
            file='reports/TF2GeneSelection.txt',sep="\t",quote=F,row.names=F)


TFs <- calcCoGAPSStat(Amean = CoGAPS_SDFIVE_3$Amean, 
                      Asd=CoGAPS_SDFIVE_3$Asd, 
                      GStoGenes=TF2Gene)
TFsUp <- apply(TFs$GSUpreg < 0.05,1,
               function(x){names(which(x))})

# perform enrichment analysis for factors regulating EGFR
EGFRTFs <- TF2Gene[grep('isoform',reverseSplit(TF2Gene)[['EGFR']],
                        value=T,invert=T)]

EGFRTFs <- sapply(EGFRTFs, intersect, row.names(CoGAPS_SDFIVE_3$Amean))
EGFRTFs <- sapply(EGFRTFs, setdiff, 'EGFR')
EGFRTFs <- EGFRTFs[sapply(EGFRTFs,length)>=5]
EGFRTFs <- EGFRTFs[sapply(EGFRTFs,length)<100]

EGFRFit <- CoGAPS_SDFIVE_3$Amean['EGFR',]%*%CoGAPS_SDFIVE_3$Pmean

apply(CoGAPS_SDFIVE_3$Amean/CoGAPS_SDFIVE_3$Asd,1,function(x){cor(x,as.numeric(CoGAPS_SDFIVE_3$Amean['EGFR',]/CoGAPS_SDFIVE_3$Asd['EGFR',]))}) -> y
y <- y[names(y)!='EGFR']

sapply(EGFRTFs,function(x){wilcoxGST(setdiff(x,'EGFR'),y)})

# plot changes in AP-2 TF targets after treatment for additional validation
delEGFR <- rep(0, length(trtSamp))
names(delEGFR) <- trtSamp

delAll <- matrix(0, nrow=nrow(HaCaT.TFGene.TrtMean),
                 ncol=length(trtSamp),
                 dimnames=list(row.names(HaCaT.TFGene.TrtMean),
                               trtSamp))
delAll <- delAll[setdiff(row.names(delAll),'EGFR'),]

for (t in c('Cetuximab','Afatinib','Gefitinib')) {
  for (c in unique(HaCaT.sampleAnnot[trtSamp,'CellLine'])) {
    samp <- trtSamp[HaCaT.sampleAnnot[trtSamp,'Treatment'] == t &
                      HaCaT.sampleAnnot[trtSamp,'CellLine']==c]
    delEGFR[samp] <- HaCaT.TFGene.TrtMean['EGFR',paste(c,'HaCaT',t,'NA')] -
      HaCaT.TFGene.TrtMean['EGFR',paste(c,'HaCaT Control NA')]
    delAll[,samp] <- HaCaT.Gene.FRMA.pSVA[row.names(delAll),samp] -
      HaCaT.TFGene.TrtMean[row.names(delAll),paste(c,'HaCaT Control NA')]
  }
}

pdf('graphs/ProbeSelectionEGFR_AP2.pdf')

probesConfirmAlpha <- sort(y[EGFRTFs[['AP-2alpha']]],decreasing=T)[1:4]

delDataPlotAlpha <- t(apply(rbind(delAll[EGFRTFs[['AP-2alpha']],],
                             EGFR=delEGFR),1,tapply,
                       HaCaT.sampleAnnot[colnames(delAll),
                                         'ExperimentalConditions'],mean))

delDataPlotAlpha <- delDataPlotAlpha[,grep('Control',colnames(delDataPlotAlpha),
                                           invert=T,value=T)]

heatmap.2(delDataPlotAlpha,
          Colv=T,
          scale='row',trace='none',col=greenred,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          labCol=sub('HaCaT-','',strspliti(colnames(delDataPlotAlpha),
                                           split=" ",i=1)),
          xlab='HaCaT Construct',
          ColSideColors=getCol(strspliti(colnames(delDataPlotAlpha), 
                                         split=" ",i=3)),
          main='AP-2alpha')



dev.off()




