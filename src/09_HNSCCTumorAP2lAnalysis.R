## load in cached data from Project Template
library('ProjectTemplate')
load.project()

## load in fRMA normalized cetuximab treated data in HNSCC tumors: 
## http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4741452/
load('~/Dropbox/SchmitzCetuximab/HumanCtxSchmidtz_09May2016.Rda')

## Finding p16 expression as a surrogate for HPV status
CDKN2AExprs <- matrix(NA, nrow=length(levels(pData(frmaData)$condition)),
                      ncol=length(unique(pData(frmaData)$Patient.number)),
                      dimnames=list(levels(pData(frmaData)$condition),
                                    unique(pData(frmaData)$Patient.number)))

for (c in row.names(CDKN2AExprs)) {
  samp <-  colnames(frmaData)[pData(frmaData)$condition == c]
  CDKN2AExprs[c,pData(frmaData)[samp,'Patient.number']] <- exprs(frmaData)[geneProbeSelect['CDKN2A'],samp]
}




pdf('graphs/SchmitzCetuximab_P16Expression.pdf',height=3.5)
barplot(CDKN2AExprs[,order(apply(CDKN2AExprs,2,mean))],beside=T, ylab=sprintf('CDKN2A mRNA expression (%s)',
                                           geneProbeSelect['CDKN2A']))
abline(h=median(apply(CDKN2AExprs,2,mean)), lty=2)
dev.off()

## make HPV-calls from p16 expression
HPV <- rep('HPVPos',ncol(CDKN2AExprs))
names(HPV) <- colnames(CDKN2AExprs)

HPVNeg <- names(which(apply(CDKN2AExprs,2,max) < median(apply(CDKN2AExprs,2,mean))))
HPV[HPVNeg] <- 'HPVNeg'
names(HPV) <- paste0('P',names(HPV))
HPVNeg <- paste0('P',HPVNeg)

## computing change in gene expression data after treatment in each patient
delExprs <- matrix(NA, nrow=nrow(frmaData) ,ncol=length(unique(pData(frmaData)$Patient.number)),
                   dimnames=list(row.names(frmaData),unique(pData(frmaData)$Patient.number)))
delExprs <- delExprs[geneProbeSelect,]
row.names(delExprs) <- names(geneProbeSelect)

pretrt <- colnames(frmaData)[pData(frmaData)$condition == 'Baseline']
trt <- colnames(frmaData)[pData(frmaData)$condition == 'Aftertreatment']

if (all(pData(frmaData)[pretrt,'Patient.number'] == pData(frmaData)[trt,'Patient.number'])) {
  message('Pre and post treatment samples are aligned')
} else {
  stop('Pre and post treatment samples are not aligned')
}

delExprs <- exprs(frmaData)[,trt] - exprs(frmaData)[,pretrt]
colnames(delExprs) <- paste0('P',pData(frmaData)[pretrt,'Patient.number'])
delExprs <- delExprs[geneProbeSelect,]
row.names(delExprs) <- names(geneProbeSelect)

## change in EGFR vs p16
pdf('graphs/EGFRvP16Exprs.pdf')
boxplot(delExprs['EGFR',HPVNeg],delExprs['EGFR',setdiff(colnames(delExprs),HPVNeg)], 
        names=paste(c('<','>'),round(median(apply(CDKN2AExprs,2,mean)),1)), 
        xlab=sprintf('p16 mRNA expression (%s)',geneProbeSelect['CDKN2A']),
        ylab=sprintf('EGFR mRNA expression post - pre-treatment (%s)',geneProbeSelect['EGFR']))
dev.off()

## load in AP-2alpha gene expression signature 
AP2Up <- read.table('reports/LINCS2000/AP2Up.grp',header=T)[,1]
AP2Down <- read.table('reports/LINCS2000/AP2Down.grp',header=T)[,1]

pdf('graphs/HumanHNSCCHeatmapAP2.pdf')
heatmap.2(delExprs[c(AP2Up,AP2Down),],
          scale='row',trace='none',col=greenred,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          labRow=c(AP2Up,AP2Down),colRow = rep(c('red','blue'),c(length(AP2Up),length(AP2Down))),
          ColSideColors = ifelse(colnames(delExprs)%in%HPVNeg,'red','blue'))
dev.off()

APStats <- matrix(NA,nrow=2,ncol=2,
                 dimnames = list(unique(HPV),c('up','down')))

pdf('graphs/HumanHNSCCTStats.pdf')
par(mfrow=c(2,1))
for (h in row.names(APStats)) {

  delAllHNSCC <- delExprs[setdiff(row.names(delExprs),'EGFR'),names(which(HPV==h))]
  delEGFRHNSCC <- delExprs['EGFR',names(which(HPV==h))]

  mm <- model.matrix(~delEGFRHNSCC)
  delFit <- eBayes(lmFit(delAllHNSCC[,names(delEGFRHNSCC)],mm))

  APStats[h,'up'] <- geneSetTest(statistics = delFit$t[,2],
                              index=intersect(AP2Up,
                                              row.names(delAllHNSCC)),
                              alternative='greater')
  APStats[h,'down'] <- geneSetTest(statistics = delFit$t[,2],
                                index=intersect(AP2Down,
                                                row.names(delAllHNSCC)),
                                alternative='less')
  
  
  plot(sort(delFit$t[,2]),type='l',ylab='t EGFR association',ylim=c(-10,10))
  
  for (p in setdiff(AP2Up,'EGFR')) {
    points(x=which(names(sort(delFit$t[,2]))  ==p),
           y=delFit$t[p,2],
           pch=19,col='red')
    text(x=which(names(sort(delFit$t[,2]))  ==p),
         y=delFit$t[p,2]+0.5,
         labels=p,col='red',srt = 90,pos=4,offset=0)
  }
  
  for (p in setdiff(AP2Down,'EGFR')) {
    points(x=which(names(sort(delFit$t[,2]))  ==p),
           y=delFit$t[p,2],
           pch=19,col='blue')
    text(x=which(names(sort(delFit$t[,2]))  ==p),
         y=delFit$t[p,2]-0.5,
         labels=p,col='blue',srt = 90,pos=2,offset=0)
  }
  
  title(sprintf('%s Up %0.2f, Down %0.2f',h,
                APStats[h,'up'],APStats[h,'down']))
  
}
dev.off()
