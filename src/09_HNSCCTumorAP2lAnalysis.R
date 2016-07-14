## load in cached data from Project Template
library('ProjectTemplate')
load.project()

## load in fRMA normalized cetuximab treated data in HNSCC tumors: 
## http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4741452/
load('HumanCtxSchmidtz_09May2016.Rda')

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



## analysis of premalignancy data

## change in gene expression for premalignancies 
## computing average of multiple pre/post treatment samples for the same patient
Premalignancy.del <- matrix(NA, nrow=length(geneProbeSelect), ncol=2,
                            dimnames = list(names(geneProbeSelect),c('F','H')))

for (p in colnames(Premalignancy.del)) {
  samp <- row.names(Premalignancy.sampleAnnot)[Premalignancy.sampleAnnot$Patient == p] 
  samp <- samp[!is.na(samp)]
  
  pre <- samp[Premalignancy.sampleAnnot[samp,'Treatment'] == 'Pre']
  post <- samp[Premalignancy.sampleAnnot[samp,'Treatment'] == 'Post']
  
  if (length(post) > 1) {
    Premalignancy.del[,p] <- apply(Premalignancy.FRMA[geneProbeSelect,post],1,mean) 
  } else {
    Premalignancy.del[,p] <- Premalignancy.FRMA[geneProbeSelect,post]
  }
  
  if (length(pre) > 1) {
    Premalignancy.del[,p] <- Premalignancy.del[,p] - apply(Premalignancy.FRMA[geneProbeSelect,pre],1,mean) 
  } else {
    Premalignancy.del[,p] <- Premalignancy.del[,p] - Premalignancy.FRMA[geneProbeSelect,pre]
  }
}
  
# change in expression between samples, ranked by EGFR expression
delDelPremalignancy <- Premalignancy.del[,names(which.max(Premalignancy.del['EGFR',]))] - 
  Premalignancy.del[,names(which.min(Premalignancy.del['EGFR',]))]


APStatsPreUp <- geneSetTest(statistics = delDelPremalignancy,
                               index=intersect(AP2Up,
                                               names(delDelPremalignancy)),
                               alternative='greater')
APStatsPreDown <- geneSetTest(statistics = delDelPremalignancy,
                                 index=intersect(AP2Down,
                                                 names(delDelPremalignancy)),
                                 alternative='less')

pdf('graphs/delDelPremalignancyAP2.pdf',height=3.5)
par(mfrow=c(1,1))
plot(sort(delDelPremalignancy),type='l',
     ylim=c(-8,8),ylab='del del mRNA expression')

for (p in AP2Up) {
  points(x=which(names(sort(delDelPremalignancy))  ==p),
         y=delDelPremalignancy[p],
         pch=19,col='red')
  text(x=which(names(sort(delDelPremalignancy))  ==p),
         y=delDelPremalignancy[p]+0.5,
         labels=p,col='red',srt = 90,pos=4,offset=0)
}

for (p in setdiff(AP2Down,'EGFR')) {
  points(x=which(names(sort(delDelPremalignancy))  ==p),
         y=delDelPremalignancy[p],
         pch=19,col='blue')
  text(x=which(names(sort(delDelPremalignancy))  ==p),
       y=delDelPremalignancy[p]-0.5,
       labels=p,col='blue',srt = 90,pos=2,offset=0)
}

title(sprintf('Up %0.2f, Down %0.2f',APStatsPreUp,APStatsPreDown))
dev.off()
