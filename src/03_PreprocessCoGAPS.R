library('ProjectTemplate')
load.project()

# compute the mean expression in each condition
HaCaT.Gene.FRMA.pSVA.Experiment <- t(apply(HaCaT.Gene.FRMA.pSVA,1,
  function(x){
    tapply(x,HaCaT.sampleAnnot[colnames(HaCaT.Gene.FRMA.pSVA),
                               'ExperimentalConditions'],mean)}))

# assign variability to be at least 5% of the signal in the data
HaCaT.Gene.FRMA.pSVA.Experiment.sd <- t(apply(HaCaT.Gene.FRMA.pSVA,1,
  function(x){
    tapply(x,HaCaT.sampleAnnot[colnames(HaCaT.Gene.FRMA.pSVA),
                               'ExperimentalConditions'],sd)}))


HaCaT.Gene.FRMA.pSVA.Experiment.sd <- pmax(HaCaT.Gene.FRMA.pSVA.Experiment.sd,
                                           0.05*HaCaT.Gene.FRMA.pSVA.Experiment)

# get the set of transcription factor targets
load('data/TRANSFAC_Genes_2014.Rda')
TFGenes <- unique(unlist(TF2Gene))
TFGenes <- intersect(TFGenes,
                     row.names(HaCaT.Gene.FRMA.pSVA.Experiment))

# select genes with enough variability
TFGenesVar <- 
  TFGenes[which(apply(HaCaT.Gene.FRMA.pSVA.Experiment.sd[TFGenes,] /
                        HaCaT.Gene.FRMA.pSVA.Experiment[TFGenes,],
                      1,max,na.rm=T)  > 0.01)]

# limit to probes with at least 0.5 FC in any conditions
fcGenes <- names(which(apply(HaCaT.Gene.FRMA.pSVA,1,
                             function(x){max(x)-min(x)}) > 0.5))
TFGenesVar <- intersect(TFGenesVar,fcGenes)

# select genes annotated to TFs with enough variability
HaCaT.TFGene.ExpMean <- HaCaT.Gene.FRMA.pSVA.Experiment[TFGenesVar,]
HaCaT.TFGene.ExpSD <- HaCaT.Gene.FRMA.pSVA.Experiment.sd[TFGenesVar,]


# get those samples that are treated only
trtExp <- colnames(HaCaT.TFGene.ExpMean)[
  strspliti(colnames(HaCaT.TFGene.ExpMean),split=" ",i=3) != "NA" &
  strspliti(colnames(HaCaT.TFGene.ExpMean),split=" ",i=4) == "NA"]

# subset data to treated samples only
HaCaT.TFGene.TrtMean <- HaCaT.TFGene.ExpMean[,trtExp]
HaCaT.TFGene.TrtSD <- HaCaT.TFGene.ExpSD[,trtExp]

ProjectTemplate::cache('HaCaT.TFGene.TrtMean')
ProjectTemplate::cache('HaCaT.TFGene.TrtSD')
