library('ProjectTemplate')
load.project()

# Perform CoGAPS analysis and establish the appropriate number of 
# patterns (unsupervised)
for (nP in 2:6) {
  
  # check whether CoGAPS has been run previously or not
  if (length(ls(pattern=sprintf('CoGAPS_SDFIVE_%d', nP))) > 0) {
    message('Using cached CoGAPS results for ', nP, ' patterns')
    next
  }
  
  message('CoGAPS for ',nP)
  assign(sprintf('CoGAPS_SDFIVE_%d',nP),
         gapsRun(D=HaCaT.TFGene.TrtMean,
                 S=HaCaT.TFGene.TrtSD,
                 nFactor=as.character(nP),
                 nEquil="50000",
                 nSample="50000"))
  ProjectTemplate::cache(sprintf('CoGAPS_SDFIVE_%d',nP))
}


pdf(sprintf('graphs/CoGAPS_SDFIVE_MeanChi2.pdf'))
plot(2:6,
     sapply(sprintf('CoGAPS_SDFIVE_%d', 2:6),
            function(x){get(x)$meanChi2}),
     xlab='Number of patterns', ylab='mean chi^2', 
     ylim=c(0,10000))
dev.off()