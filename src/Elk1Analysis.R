library('ProjectTemplate')
load.project()

load('data/TRANSFAC_Genes_2014.Rda')

TF2Gene <- sapply(TF2Gene,intersect,
                  row.names(CoGAPS_SDFIVE_3$Amean))

set.seed(1234)
Elk1Stats <- calcCoGAPSStat(CoGAPS_SDFIVE_3$Amean,CoGAPS_SDFIVE_3$Asd,TF2Gene['Elk-1'])

pdf('graphs/Elk1Targets_Pattern2.pdf',width=3.5,height=3.5)
for (i in 2) {
  plot(sort(CoGAPS_SDFIVE_3$Amean[,i]),type='l')
  polygon(c(1:nrow(CoGAPS_SDFIVE_3$Amean),
            seq(from=nrow(CoGAPS_SDFIVE_3$Amean),to=1)),
          c(sort(CoGAPS_SDFIVE_3$Amean[,i]) + 
              CoGAPS_SDFIVE_3$Asd[names(sort(CoGAPS_SDFIVE_3$Amean[,i])),3],
            sort(CoGAPS_SDFIVE_3$Amean[,i],decreasing=T) - 
              CoGAPS_SDFIVE_3$Asd[names(sort(CoGAPS_SDFIVE_3$Amean[,i],decreasing=T)),3]),col='grey',
          border='grey')
  lines(sort(CoGAPS_SDFIVE_3$Amean[,i]),type='l',lwd=2)
  for (g in TF2Gene[['Elk-1']]) {
    EGFRidx <- which(names(sort(CoGAPS_SDFIVE_3$Amean[,i]))==g)
    errbar(EGFRidx,CoGAPS_SDFIVE_3$Amean[g,i],errbar.col = 'red',
           CoGAPS_SDFIVE_3$Amean[g,i]+CoGAPS_SDFIVE_3$Asd[g,i],
           CoGAPS_SDFIVE_3$Amean[g,i]-CoGAPS_SDFIVE_3$Asd[g,i],
           col='red',pch='19',add=T,lwd=2)
  }
  title(sprintf('Association of Elk-1 Targets with Pattern 2\n p=%0.2f',
                Elk1Stats$GSUpreg[2,]))
}
dev.off()