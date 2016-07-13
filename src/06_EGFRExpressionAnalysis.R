library('ProjectTemplate')
load.project()

### association of EGFR with each of the CoGAPS patterns
g <- 'EGFR'

# compute simplicty statistics for CoGAPS Patterns


pdf('graphs/CoGAPS_RankEGFR.pdf', height=3.5)
par(mfrow=c(1,3))
for (i in 1:3) {
  plot(sort(CoGAPS_SDFIVE_3$Amean[,i]),type='l',
       xlab='Gene rank',
       ylab=paste('Amplitude for pattern',i))
  polygon(c(1:nrow(CoGAPS_SDFIVE_3$Amean),
            seq(from=nrow(CoGAPS_SDFIVE_3$Amean),to=1)),
          c(sort(CoGAPS_SDFIVE_3$Amean[,i]) + 
              CoGAPS_SDFIVE_3$Asd[names(sort(CoGAPS_SDFIVE_3$Amean[,i])),3],
            sort(CoGAPS_SDFIVE_3$Amean[,i],decreasing=T) - 
              CoGAPS_SDFIVE_3$Asd[names(sort(CoGAPS_SDFIVE_3$Amean[,i],
                                             decreasing=T)),3]),
          col='grey', border='grey')
  lines(sort(CoGAPS_SDFIVE_3$Amean[,i]),type='l',lwd=2)
  EGFRidx <- which(names(sort(CoGAPS_SDFIVE_3$Amean[,i]))==g)
  errbar(EGFRidx,CoGAPS_SDFIVE_3$Amean[g,i],errbar.col = 'red',
         CoGAPS_SDFIVE_3$Amean[g,i]+CoGAPS_SDFIVE_3$Asd[g,i],
         CoGAPS_SDFIVE_3$Amean[g,i]-CoGAPS_SDFIVE_3$Asd[g,i],
         col='red',pch='19',add=T,lwd=2,cap=0)
}
dev.off()

# differential expression analysis for treatment overall
trtSamp <- row.names(HaCaT.sampleAnnot)[
  paste(HaCaT.sampleAnnot$Treatment,HaCaT.sampleAnnot$siRNA) %in% 
    c('Control NA', 'Cetuximab NA', 'Afatinib NA', 'Gefitinib NA')]
trtSamp <- intersect(trtSamp, colnames(HaCaT.Gene.FRMA.pSVA))

trtModel <- with(HaCaT.sampleAnnot[trtSamp,], 
                 model.matrix(~0+Treatment+CellLine))
colnames(trtModel) <- gsub('-','',colnames(trtModel))
trtLm <- lmFit(HaCaT.Gene.FRMA.pSVA[,trtSamp],trtModel)
trtCon <- eBayes(contrasts.fit(trtLm,
                               makeContrasts(contrasts=c("TreatmentCetuximab-TreatmentControl",
                                                         "TreatmentGefitinib-TreatmentControl",
                                                         "TreatmentAfatinib-TreatmentControl"),
                                             levels=trtModel)))


# compute the average EGFR family expression
  
EGFRExprs <- matrix(NA, nrow=4, ncol=6,
                    dimnames=list(c('HaCaT-Mock','HaCaT-EGFR',
                                    'HaCaT-HRAS','HaCaT-PIK3CA'),
                                  c('Control NA', 'Cetuximab NA', 
                                    'Gefitinib NA', 'Afatinib NA',
                                    'NA Scramble', 'NA EGFR')))
EGFRMin <- EGFRMax <- EGFRP <- EGFRExprs

for (c in row.names(EGFRExprs)) {
  
  trtSamp <- row.names(HaCaT.sampleAnnot)[
    paste(HaCaT.sampleAnnot$Treatment,HaCaT.sampleAnnot$siRNA) %in% 
      c('Control NA', 'Cetuximab NA', 'Afatinib NA', 'Gefitinib NA', 'NA Scramble', 'NA EGFR')]
  trtSamp <- trtSamp[HaCaT.sampleAnnot[trtSamp,'CellLine']==c]
  trtSamp <- intersect(trtSamp, colnames(HaCaT.Gene.FRMA.pSVA))
  
  Treatment <- paste0(HaCaT.sampleAnnot[trtSamp,'Treatment'],
                     HaCaT.sampleAnnot[trtSamp,'siRNA'])
  
  trtModel <- model.matrix(~0+ Treatment)
  colnames(trtModel) <- gsub('-','',colnames(trtModel))
  trtLm <- lmFit(HaCaT.Gene.FRMA.pSVA[,trtSamp],trtModel)
  trtCon <- eBayes(contrasts.fit(trtLm,
                                 makeContrasts(contrasts=c("TreatmentCetuximabNA-TreatmentControlNA",
                                                           "TreatmentGefitinibNA-TreatmentControlNA",
                                                           "TreatmentAfatinibNA-TreatmentControlNA",
                                                           "TreatmentNAEGFR - TreatmentNAScramble"),
                                               levels=trtModel)))
  
  
  for (e in colnames(EGFRExprs)) {
    samp <- colnames(HaCaT.Gene.FRMA.pSVA)
    samp <- samp[HaCaT.sampleAnnot[samp,'CellLine']==c]
    samp <- samp[paste(HaCaT.sampleAnnot[samp,'Treatment'],
                       HaCaT.sampleAnnot[samp,'siRNA'])==e]
    EGFRExprs[c,e] <- mean(HaCaT.Gene.FRMA.pSVA[g,samp])
    EGFRMin[c,e] <- min(HaCaT.Gene.FRMA.pSVA[g,samp])
    EGFRMax[c,e] <- max(HaCaT.Gene.FRMA.pSVA[g,samp])
    
    if (e!= 'Control NA' & e!= 'NA Scramble') {
      EGFRP[c,e] <- trtCon['EGFR',grep(sub(" ","",e),colnames(trtCon),value=T)]$p.value
    }
  }
}

# univariate p-values by cell line and treatment relative to control
EGFRP <- EGFRP[,!apply(is.na(EGFRP),2,any)]
pdf('graphs/EGFRPValues.pdf')
grid.table(signif(EGFRP,2))
dev.off()


EGFRStats <- rep(NA, 4)
names(EGFRStats) <- c('Cetuximab NA vs Control NA', 
                      'Gefitinib NA vs Control NA',
                      'Afatinib NA vs Control NA',
                      'NA EGFR vs NA Scramble')
for (t in names(EGFRStats)) {
  grps <- strsplit(t,split=" vs ")[[1]]
  
  samp <- colnames(HaCaT.Gene.FRMA.pSVA)
  samp <- samp[paste(HaCaT.sampleAnnot[samp,'Treatment'],
                     HaCaT.sampleAnnot[samp,'siRNA']) %in% grps]
  
  EGFRStats[t] <- 
    anova(lm(HaCaT.Gene.FRMA.pSVA[g,samp]~
               HaCaT.sampleAnnot[samp,'CellLine'] + 
               paste(HaCaT.sampleAnnot[samp,'Treatment'],
                     HaCaT.sampleAnnot[samp,'siRNA'])),
          lm(HaCaT.Gene.FRMA.pSVA[g,samp]~
               HaCaT.sampleAnnot[samp,'CellLine']))$`Pr(>F)`[2]
}


pdf(sprintf('graphs/%sExpressionTreatment.pdf',g),height=3.5)
par(mfrow=c(1,3))

plot(c(0,20),c(5,12),col='white',
     ylab=sprintf('EGFR mRNA expression (%s)',geneProbeSelect[g]),
     xlab='',axes=F)
for (i in 1:4) {
  errbar(seq(from=4*(i-1) + .5 + i,to=4*i + .5 + (i-1)),EGFRExprs[i,1:4],
         EGFRMax[i,1:4],EGFRMin[i,1:4],add=T,
         pch=getPCH(row.names(EGFRExprs)[i]), lwd=2,
         col=getCol(strspliti(colnames(EGFRExprs)[1:4],split=" ",i=1)),
         errbar.col=getCol(strspliti(colnames(EGFRExprs)[1:4],split=" ",i=1)))
}
title(sprintf('CTX v Control %0.1e\nGFT v Control %0.1e\nAFT v Control %0.1e',
              topTable(trtCon,number=Inf,
                       coef="TreatmentCetuximab-TreatmentControl")[g,'adj.P.Val'], 
              topTable(trtCon,number=Inf,
                       coef="TreatmentGefitinib-TreatmentControl")[g,'adj.P.Val'],
              topTable(trtCon,number=Inf,
                       coef="TreatmentAfatinib-TreatmentControl")[g,'adj.P.Val']))
axis(2,las=1)
axis(1,at=sapply(1:4,function(i){mean(seq(from=4*(i-1) + .5 + i,
                                          to=4*i + .5 + (i-1)))}),
     labels=row.names(EGFRExprs),
     las=2)

plot(c(0,12),c(5,8),col='white',
     ylab=sprintf('EGFR mRNA expression (%s)',geneProbeSelect[g]),
     xlab='',axes=F)
for (i in 1:4) {
  errbar(seq(from=2*(i-1)+.5+i,to=2*i+.5+(i-1)),
         EGFRExprs[i,5:6], EGFRMax[i,5:6], EGFRMin[i,5:6], 
         add=T, pch=getPCH(row.names(EGFRExprs)[i]), 
         lwd=2,col=c(grey(0.25),'orange'),
         errbar.col=c(grey(0.25),'orange'))
}
axis(2,las=1)
axis(1,at=sapply(1:4,function(i){mean(seq(from=2*(i-1)+.5+i,
                                          to=2*i+.5+(i-1)))}),
     labels=row.names(EGFRExprs),
     las=2)
title(sprintf('siEGFR v scramble %0.1e', 
              EGFRStats['NA EGFR vs NA Scramble']))

#dev.off()

delEGFR <- sweep(EGFRExprs,1,EGFRExprs[,'Control NA'])
#pdf(sprintf('graphs/%svsDrugSens.pdf',g))
plot(100*DrugSens[,'Cetuximab'], 
     delEGFR[row.names(DrugSens),'Cetuximab NA'], 
     col=getCol('Cetuximab'), pch=getPCH(row.names(EGFRExprs)),
     ylim=c(0,max(delEGFR)),xlim=c(0,100),
     ylab=sprintf('%s exprs in Trt - %s exprs in Control',g,g),
     xlab='Survival (%)')
points(100*DrugSens[,'Gefitinib'], 
       delEGFR[row.names(DrugSens),'Gefitinib NA'], 
       col=getCol('Gefitinib'), pch=getPCH(row.names(EGFRExprs)))
points(100*DrugSens[,'Afatinib'], 
       delEGFR[row.names(DrugSens),'Afatinib NA'], 
       col=getCol('Afatinib'), pch=getPCH(row.names(EGFRExprs)))
title(sprintf('Cor %0.1f p %0.2f',
              cor(as.numeric(DrugSens),
                  as.numeric(delEGFR[row.names(DrugSens),
                                     paste(colnames(DrugSens),'NA')])),
              cor.test(as.numeric(DrugSens),
                       as.numeric(delEGFR[row.names(DrugSens),
                                          paste(colnames(DrugSens),'NA')]))$p.value))
abline(lm(as.numeric(delEGFR[row.names(DrugSens),
                             paste(colnames(DrugSens),'NA')])
          ~ as.numeric(100*DrugSens)),
       lty=2, lwd=2, col='black')
dev.off()

pdf('graphs/Legend.pdf')
plot(c(0,1),c(0,1),axes=F,xlab="",ylab="",col='white')
legend('topright',pch='-',
       col=c(getCol(c('Control','Cetuximab', 'Gefitinib','Afatinib')),
                    grey(0.25),'orange'),
       legend=c('Control','Cetuximab', 'Gefitinib','Afatinib',
                'siScramble','siEGFR'))
legend('topleft',
       pch=c(1,0,5,2),
       legend=row.names(EGFRExprs))
dev.off()
  



