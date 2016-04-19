### fRMA normalize the HaCaT data
HaCaT.FRMA.BATCH1 <- frma(HaCaT.CEL.BATCH1)
HaCaT.FRMA.BATCH2 <- frma(HaCaT.CEL.BATCH2)
HaCaT.FRMA.BATCH3 <- frma(HaCaT.CEL.BATCH3)

### Merge data from two batches
HaCaT.FRMA <- cbind(exprs(HaCaT.FRMA.BATCH1),
                    exprs(HaCaT.FRMA.BATCH2),
                    exprs(HaCaT.FRMA.BATCH3))

# save frma data in the cache
ProjectTemplate::cache('HaCaT.FRMA')

### Get the sample annotations
sampleAnnot <- data.frame(SampleName=strspliti(colnames(HaCaT.FRMA),
                                               split="-",3),
                          stringsAsFactors=F)
row.names(sampleAnnot) <- colnames(HaCaT.FRMA)

# extract the overexpression construct of the cell line
sampleAnnot$CellLine <- rep('HaCaT-Mock',nrow(sampleAnnot))
sampleAnnot$CellLine[grep('SCC[12]',sampleAnnot$SampleName)] <- 
  grep('SCC[12]',sampleAnnot$SampleName, value=T)
sampleAnnot$CellLine[sampleAnnot$SampleName=='FADU'] <- 'FaDu'
sampleAnnot$CellLine[sampleAnnot$SampleName=='HACAT'] <- 'HaCaT-Parent'

sampleAnnot$CellLine[grep('EG',sampleAnnot$SampleName)] <- 'HaCaT-EGFR'
sampleAnnot$CellLine[grep('HR',sampleAnnot$SampleName)] <- 'HaCaT-HRAS'
sampleAnnot$CellLine[grep('PI',sampleAnnot$SampleName)] <- 'HaCaT-PIK3CA'

# get the cell line name
sampleAnnot$CellType <- strspliti(sampleAnnot$CellLine,split="-",1)

# get the treatment type
sampleAnnot$Treatment <- rep(NA, nrow(sampleAnnot))
sampleAnnot$Treatment[grep('A[FP]',sampleAnnot$SampleName)] <- 'Afatinib'
sampleAnnot$Treatment[grep('ATF',sampleAnnot$SampleName)] <- 'Afatinib'
sampleAnnot$Treatment[substr(sampleAnnot$SampleName,5,6)=='GF'] <- 'Gefitinib'
sampleAnnot$Treatment[substr(sampleAnnot$SampleName,7,8)=='GF'] <- 'Gefitinib'
sampleAnnot$Treatment[substr(sapply(strsplit(sampleAnnot$SampleName,split="_"),
       function(x){x[length(x)]}),1,2)=="GF"] <- 'Gefitinib'
sampleAnnot$Treatment[grep('GTF',sampleAnnot$SampleName)] <- 'Gefitinib'

sampleAnnot$Treatment[grep('C[TX]',sampleAnnot$SampleName)] <- 'Cetuximab'
sampleAnnot$Treatment[grep('SC_?CO',sampleAnnot$SampleName)] <- 'Control'
sampleAnnot$Treatment[grep('[Cc]t[rn]l',sampleAnnot$SampleName)] <- 'Control'

# determine whether there was siRNA knock down or scramble
sampleAnnot$siRNA <- rep(NA, nrow(sampleAnnot))
sampleAnnot$siRNA[grep('KD',sampleAnnot$SampleName)] <- 'EGFR'
sampleAnnot$siRNA[grep('SC[0-9]',sampleAnnot$SampleName)] <- 'Scramble'
sampleAnnot$siRNA[grep('_SC',sampleAnnot$SampleName)] <- 'Scramble'
sampleAnnot$siRNA[substr(sampleAnnot$SampleName,5,6)=='SC'] <- 'Scramble'

# determine the amount of time at which the cell lines were processed
sampleAnnot$Timing <- rep(24, nrow(sampleAnnot))
sampleAnnot$Timing[!is.na(sampleAnnot$Treatment)] <- 29
sampleAnnot$Timing[sampleAnnot$Treatment=='Control'] <- 29

# add the cell line batch
sampleAnnot$Batch <- substr(row.names(sampleAnnot),1,1)

# standardized notation for the experimental conditions of the cell lines in 
sampleAnnot$ExperimentalConditions <- paste(sampleAnnot$CellLine,
                                            sampleAnnot$CellType,
                                            sampleAnnot$Treatment,
                                            sampleAnnot$siRNA)

# save data for caching
HaCaT.sampleAnnot <- sampleAnnot
ProjectTemplate::cache('HaCaT.sampleAnnot')

# drug sensitivity data
for (f in list.files('data',pattern='ColonyFormation')) {
  d <- strspliti(f,split="_",i=2)
  
  CF <- read.table(file.path('data',f),
                   header=T,sep=",")
  
  
  if (any(colnames(CF)=='cell.line')) {
    CF$Experiment.Title <- CF$cell.line
  }
  
  if (any(colnames(CF)=='total.area')) {
    CF$Total.area <- CF$total.area
  }
  
  CF$Treatment.nM <- as.numeric(sub('PBS','0',
                                    sub('DMSO','0',sub('nM','',CF$Treatment))))
  
  TrtMean <- TrtMin <- TrtMax <- 
    TotalAreaMean <- TotalAreaMin <- TotalAreaMax <- 
    matrix(NA, nrow=length(unique(CF$Experiment.Title)),
           ncol=length(unique(CF$Treatment.nM)),
           dimnames=list(unique(CF$Experiment.Title),
                         unique(CF$Treatment.nM)))
  
  for (c in unique(CF$Experiment.Title)) {
    TrtMean[c,] <- 
      tapply(CF$TA.prop[CF$Experiment.Title==c],
             CF$Treatment.nM[CF$Experiment.Title==c],
             mean)
    
    TrtMin[c,] <- 
      tapply(CF$TA.prop[CF$Experiment.Title==c],
             CF$Treatment.nM[CF$Experiment.Title==c],
             min)
    
    TrtMax[c,] <- 
      tapply(CF$TA.prop[CF$Experiment.Title==c],
             CF$Treatment.nM[CF$Experiment.Title==c],
             max)
    
    TotalAreaMean[c,] <- 
      tapply(CF$Total.area[CF$Experiment.Title==c],
             CF$Treatment.nM[CF$Experiment.Title==c],
             mean)
    
    TotalAreaMin[c,] <- 
      tapply(CF$Total.area[CF$Experiment.Title==c],
             CF$Treatment.nM[CF$Experiment.Title==c],
             min)
    
    TotalAreaMax[c,] <- 
      tapply(CF$Total.area[CF$Experiment.Title==c],
             CF$Treatment.nM[CF$Experiment.Title==c],
             max)
  }
  
  assign(sprintf('HaCaT.%s.TrtMean',d),TrtMean)
  assign(sprintf('HaCaT.%s.TrtMin',d),TrtMin)
  assign(sprintf('HaCaT.%s.TrtMax',d),TrtMax)
  
  assign(sprintf('HaCaT.%s.TotalAreaMean',d),TotalAreaMean)
  assign(sprintf('HaCaT.%s.TotalAreaMin',d),TotalAreaMin)
  assign(sprintf('HaCaT.%s.TotalAreaMax',d),TotalAreaMax)
  
  ProjectTemplate::cache(sprintf('HaCaT.%s.TrtMean',d))
  ProjectTemplate::cache(sprintf('HaCaT.%s.TrtMin',d))
  ProjectTemplate::cache(sprintf('HaCaT.%s.TrtMax',d))
  
  ProjectTemplate::cache(sprintf('HaCaT.%s.TotalAreaMean',d))
  ProjectTemplate::cache(sprintf('HaCaT.%s.TotalAreaMin',d))
  ProjectTemplate::cache(sprintf('HaCaT.%s.TotalAreaMax',d))
  
}

DrugSens <- DrugSens.Max <- DrugSens.Min <- matrix(NA_real_, nrow=4, ncol=3, 
  dimnames=list(c('HaCaT-Mock','HaCaT-EGFR','HaCaT-HRAS','HaCaT-PIK3CA'),
    c('Cetuximab','Gefitinib','Afatinib')))
DrugSens[,'Cetuximab'] <- HaCaT.CTX.TrtMean[,'100']
DrugSens[,'Gefitinib'] <- HaCaT.GFT.TrtMean[,'100']
DrugSens[,'Afatinib']  <- HaCaT.AFT.TrtMean[,'10'] 

DrugSens.Max[,'Cetuximab'] <- HaCaT.CTX.TrtMax[,'100']
DrugSens.Max[,'Gefitinib'] <- HaCaT.GFT.TrtMax[,'100']
DrugSens.Max[,'Afatinib'] <- HaCaT.AFT.TrtMax[,'10']

DrugSens.Min[,'Cetuximab'] <- HaCaT.CTX.TrtMin[,'100']
DrugSens.Min[,'Gefitinib'] <- HaCaT.GFT.TrtMin[,'100']
DrugSens.Min[,'Afatinib'] <- HaCaT.AFT.TrtMin[,'10']


ProjectTemplate::cache('DrugSens')
ProjectTemplate::cache('DrugSens.Max')
ProjectTemplate::cache('DrugSens.Min')
