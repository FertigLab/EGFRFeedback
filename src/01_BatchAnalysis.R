library('ProjectTemplate')
load.project()

# Asses whether there is an apparent batch effect in the data
pdf('graphs/AllHaCaTDataClust.pdf',width=12)
for (c in colnames(HaCaT.sampleAnnot)) {
  plotColoredClusters(standard.pearson(HaCaT.FRMA),
                      labs=paste(HaCaT.sampleAnnot[colnames(HaCaT.FRMA),
                                                   'ExperimentalConditions'],
                                 HaCaT.sampleAnnot[colnames(HaCaT.FRMA),
                                                   'Batch']),
                      cols=as.character(as.numeric(as.factor(HaCaT.sampleAnnot[colnames(HaCaT.FRMA),c]))), 
                      xlab='',cex=0.5)
  legend('topleft', fill=as.character(1:length(unique(HaCaT.sampleAnnot[,c]))), 
         title=c,
         legend=levels(factor(HaCaT.sampleAnnot[colnames(HaCaT.FRMA),c])))
}
dev.off()

# Asses whether there is an apparent batch effect in the data
pdf('graphs/AllDataClust.pdf',width=12)
plotColoredClusters(standard.pearson(cbind(HaCaT.FRMA,
                                           HNSCC.CellLines.FRMA,
                                           Premalignancy.FRMA)),
                    labs=c(HaCaT.sampleAnnot[colnames(HaCaT.FRMA),
                                           'ExperimentalConditions'],
                           HNSCC.CellLines.sampleAnnot[colnames(HNSCC.CellLines.FRMA),'CellLine_Sen'],
                           Premalignancy.sampleAnnot[colnames(Premalignancy.FRMA),'coreID']),
                    cols=c(as.character(as.numeric(as.factor(HaCaT.sampleAnnot[colnames(HaCaT.FRMA),'Batch']))),
                           rep(c('4','5'),
                               c(ncol(HNSCC.CellLines.FRMA),ncol(Premalignancy.FRMA)))))

dev.off()

# batch correct data using pSVA
HaCaT.FRMA.pSVA <- psva(dat = HaCaT.FRMA, 
                        batch = HaCaT.sampleAnnot[colnames(HaCaT.FRMA),'Batch'])

pdf('graphs/pSVAAllDataClust.pdf',width=12)
for (c in colnames(HaCaT.sampleAnnot)) {
  plotColoredClusters(standard.pearson(HaCaT.FRMA.pSVA),
                      labs=paste(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),
                                                   'ExperimentalConditions'],
                                 HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),
                                                   'Batch']),
                      cols=as.character(as.numeric(as.factor(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),c]))), 
                      xlab='',cex=0.5)
  legend('topleft', fill=as.character(1:length(unique(HaCaT.sampleAnnot[,c]))), 
         title=c,
         legend=levels(factor(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),c])))
}
dev.off()

# filter samples used purely for batch correction
HaCaT.FRMA.pSVA <- HaCaT.FRMA.pSVA[,row.names(HaCaT.sampleAnnot)[HaCaT.sampleAnnot$CellLine %in% paste('HaCaT',c('Mock','EGFR','HRAS','PIK3CA'),sep="-")]]

HaCaT.diff <- matrix(NA, nrow=nrow(HaCaT.FRMA.pSVA), ncol=ncol(HaCaT.FRMA.pSVA),
                     dimnames=dimnames(HaCaT.FRMA.pSVA))
for (s in colnames(HaCaT.diff)) {
  e <- HaCaT.sampleAnnot[s,'ExperimentalConditions']
  sother <- setdiff(row.names(HaCaT.sampleAnnot)[which(HaCaT.sampleAnnot$ExperimentalConditions==e)],s)
  if (length(sother)<1) next

  if (length(sother) > 1) {
    sdiff <- sweep(HaCaT.FRMA.pSVA[,sother],1,HaCaT.FRMA.pSVA[,s])
    HaCaT.diff[,s] <- sdiff[,names(which.min(apply(sdiff,2,function(x){sum(x^2)})))]
  } else {
    HaCaT.diff[,s] <- HaCaT.FRMA.pSVA[,sother] - HaCaT.FRMA.pSVA[,s]
  }
}

# remove those samples which are outliers relative to their replicates
HaCaT.euclid <- apply(HaCaT.diff,2,function(x){sqrt(sum(x^2))})
sampFilt <- setdiff(colnames(HaCaT.FRMA.pSVA),
                    names(which(HaCaT.euclid > sd(HaCaT.euclid,na.rm=T)*3)))
HaCaT.FRMA.pSVA <- HaCaT.FRMA.pSVA[,sampFilt]
pdf('graphs/pSVAFilterDataClust.pdf',width=12)
for (c in colnames(HaCaT.sampleAnnot)) {
 plotColoredClusters(standard.pearson(HaCaT.FRMA.pSVA),
                     labs=paste(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),
                                                  'ExperimentalConditions'],
                                HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),
                                                  'Batch']),
                     cols=as.character(as.numeric(as.factor(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),c]))), 
                     xlab='',cex=0.5)
 legend('topleft', fill=as.character(1:length(unique(HaCaT.sampleAnnot[,c]))), 
        title=c,
        legend=levels(factor(HaCaT.sampleAnnot[colnames(HaCaT.FRMA.pSVA),c])))
}
dev.off()

ProjectTemplate::cache('HaCaT.FRMA.pSVA')

