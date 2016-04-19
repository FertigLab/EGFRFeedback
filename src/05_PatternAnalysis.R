library('ProjectTemplate')
load.project()

pdf('graphs/DoseResponse_TotalArea.pdf')
par(mfrow=c(1,3))

drugs <- c('Cetuximab', 'Gefitinib', 'Afatinib')
names(drugs) <- c('CTX','GFT', 'AFT')


for (d in names(drugs)) {
  drugMean <- get(paste('HaCaT',d,'TotalAreaMean',
                        sep='.'))
  drugMin <- get(paste('HaCaT',d,'TotalAreaMin',
                       sep='.'))
  drugMax <- get(paste('HaCaT',d,'TotalAreaMax',
                       sep='.'))
  
  errbar(x=1:ncol(drugMean),
         y=drugMean['mock',],
         yplus=drugMax['mock',],
         yminus=drugMin['mock',],
         pch=getPCH('HaCaT-Mock'),type='o',axes=F,
         xlab=paste(drugs[d], 
                    'concentration (nM)'),
         ylab='Total area',
         ylim=c(0,1.1*max(drugMax)))
  errbar(x=1:ncol(drugMean),
         y=drugMean['EGFR',],
         yplus=drugMax['EGFR',],
         yminus=drugMin['EGFR',],
         add=T,pch=getPCH('HaCaT-EGFR'),type='o',col='blue',
         errbar.col = 'blue')
  errbar(x=1:ncol(drugMean),
         y=drugMean['HRAS',],
         yplus=drugMax['HRAS',],
         yminus=drugMin['HRAS',],
         add=T,pch=getPCH('HaCaT-HRAS'),type='o',col='red',
         errbar.col = 'red')
  errbar(x=1:ncol(drugMean),
         y=drugMean['PI3K',],
         yplus=drugMax['PI3K',],
         yminus=drugMin['PI3K',],
         add=T,pch=getPCH('HaCaT-PIK3CA'),type='o',col='green',
         errbar.col = 'green')
  axis(2,las=1)
  axis(1,at = 1:ncol(drugMean),
       labels = colnames(drugMean))
  legend('topright',
         col=c('black','blue','red','green'),
         pch=getPCH(c('HaCaT-Mock','HaCaT-EGFR',
                      'HaCaT-HRAS', 'HaCaT-PIK3CA')),
         legend=c('HaCaT-Mock','HaCaT-EGFR',
                  'HaCaT-HRAS', 'HaCaT-PIK3CA'))
}
dev.off()

pdf('graphs/DoseResponse.pdf')
par(mfrow=c(1,3))

drugs <- c('Cetuximab', 'Gefitinib', 'Afatinib')
names(drugs) <- c('CTX','GFT', 'AFT')


for (d in names(drugs)) {
  drugMean <- get(paste('HaCaT',d,'TrtMean',
                        sep='.'))
  drugMin <- get(paste('HaCaT',d,'TrtMin',
                       sep='.'))
  drugMax <- get(paste('HaCaT',d,'TrtMax',
                       sep='.'))
  errbar(x=1:ncol(drugMean),
         y=100*drugMean['mock',],
         yplus=100*drugMax['mock',],
         yminus=100*drugMin['mock',],
         pch=getPCH('HaCaT-Mock'),type='o',axes=F,
         xlab=paste(drugs[d], 
                    'concentration (nM)'),
         ylab='% Total area (relative to control)',
         ylim=c(0,120))
  errbar(x=1:ncol(drugMean),
         y=100*drugMean['EGFR',],
         yplus=100*drugMax['EGFR',],
         yminus=100*drugMin['EGFR',],
         add=T,pch=getPCH('HaCaT-EGFR'),type='o',col='blue',
         errbar.col = 'blue')
  errbar(x=1:ncol(drugMean),
         y=100*drugMean['HRAS',],
         yplus=100*drugMax['HRAS',],
         yminus=100*drugMin['HRAS',],
         add=T,pch=getPCH('HaCaT-HRAS'),type='o',col='red',
         errbar.col = 'red')
  errbar(x=1:ncol(drugMean),
         y=100*drugMean['PI3K',],
         yplus=100*drugMax['PI3K',],
         yminus=100*drugMin['PI3K',],
         add=T,pch=getPCH('HaCaT-PIK3CA'),type='o',col='green',
         errbar.col = 'green')
  axis(2,at=seq(from=0,to=120,by=20),
       labels=c(seq(from=0,to=100,by=20),''),las=1)
  axis(1,at = 1:ncol(drugMean),
       labels = colnames(drugMean))
  legend('topright',
         col=c('black','blue','red','green'),
         pch=getPCH(c('HaCaT-Mock','HaCaT-EGFR',
                      'HaCaT-HRAS', 'HaCaT-PIK3CA')),
         legend=c('HaCaT-Mock','HaCaT-EGFR',
                  'HaCaT-HRAS', 'HaCaT-PIK3CA'))
}
dev.off()

# plot of relative survival rates
pdf('graphs/DrugSens.pdf')
plot(c(0,4),c(0,120),col='white',axes=F,
     xlab='', ylab='Relative Cell survival rate (%)')
axis(2)
axis(1, at=seq(from=0.3,to=3.3),
     labels=row.names(DrugSens),las=2)
for (i in 1:ncol(DrugSens)) {
  errbar(seq(from=i*.2,to=3+i*.2),
         100*DrugSens[,i],
         100*DrugSens.Min[,i],
         100*DrugSens.Max[,i],
         col=getCol(colnames(DrugSens)[i]),
         errbar.col = getCol(colnames(DrugSens)[i]),
         pch=getPCH(row.names(DrugSens)),add=T,axes=F,ylab='',xlab='')
}

legend('topleft', legend=colnames(DrugSens),
       pch=19, col=getCol(colnames(DrugSens)))
dev.off()



# heatmap of all data (w/ hierarchical clustering)
pdf('graphs/HeatmapTFGene.pdf')
heatmap.2(HaCaT.TFGene.TrtMean,scale='row',
          trace='none',col=greenred,labRow=NA,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          labCol=sub('HaCaT-','',
                     strspliti(colnames(HaCaT.TFGene.TrtMean),
                               split=" ",i=1)),
          xlab='HaCaT Construct',
          ColSideColors=getCol(strspliti(colnames(HaCaT.TFGene.TrtMean), 
                                         split=" ",i=3)))
dev.off()

# plots of CoGAPS patterns
pdf('graphs/CoGAPS_SDFIVE_Patterns.pdf',width=5) #,height=3.5)

par(mfrow=c(2,2))
for (p in 1:3) {
  plot(c(0,4),
       c(0,
         max(CoGAPS_SDFIVE_3$Pmean[p,]+
               1.5*CoGAPS_SDFIVE_3$Psd[p,])),
       col='white',axes=F,
       xlab='', ylab=paste('Pattern',p))
  axis(2)
  axis(1, at=seq(from=0.3,to=3.3),
       labels=gsub('HaCaT-','',row.names(DrugSens)),las=2)
  for (i in 1:3) {
    exp <- colnames(DrugSens)[i] #c('Control',colnames(DrugSens))[i]
    points(x=seq(from=(i-1)*.3,to=3+(i-1)*.3),
           y=CoGAPS_SDFIVE_3$Pmean[p,paste(row.names(DrugSens),'HaCaT',exp,'NA')],
           col=getCol(exp),
           pch=getPCH(row.names(DrugSens)))
  }
  
  for (i in 1:nrow(DrugSens)) {
    
    segments(x0=(i-1),x1=(i-1)+.6,col="grey",
             y0=mean(CoGAPS_SDFIVE_3$Pmean[p,paste(row.names(DrugSens)[i],'HaCaT',colnames(DrugSens),'NA')]), 
             lty=2)
    

    arrows(x0=(i-1)+.3,y0 = CoGAPS_SDFIVE_3$Pmean[p,paste(row.names(DrugSens)[i],'HaCaT','Control','NA')],
           y1=mean(CoGAPS_SDFIVE_3$Pmean[p,paste(row.names(DrugSens)[i],'HaCaT',colnames(DrugSens),'NA')]),
           lty=2,col="grey",length=0.1)
    
    points(x=(i-1)+.3,y=CoGAPS_SDFIVE_3$Pmean[p,paste(row.names(DrugSens)[i],'HaCaT','Control','NA')],
           pch=getPCH(row.names(DrugSens)[i]))
  }
  
  legend('bottomleft', legend=c('Control',colnames(DrugSens),'Treated','Average change'),
         pch="-", col=c(getCol(c('Control',colnames(DrugSens))),'grey','grey'))
}

dev.off()


pdf('graphs/Pattern2_DrugSens.pdf',width=3.5,height=3.5)
CGTrtP <- CoGAPS_SDFIVE_3$Pmean[2,grep('Control',colnames(CoGAPS_SDFIVE_3$Pmean),invert=T)]

DrugSensP <- rep(NA, length(CGTrtP))
names(DrugSensP) <- names(CGTrtP)
for (i in names(DrugSensP)) {
  DrugSensP[i] <- DrugSens[strsplit(i,split=" ")[[1]][1],
                           strsplit(i,split=" ")[[1]][3]]
  CGTrtP[i] <- -CoGAPS_SDFIVE_3$Pmean[2,sub(strspliti(i,split=" ",i=3),
                                            'Control',i)] + 
    CGTrtP[i] 
}

plot(DrugSensP, CGTrtP, xlab='Cell Survival (%)',
     ylab='Pattern 2 in treatment - control',
     col=getCol(strspliti(names(DrugSensP),split=" ",i=3)),
     pch=getPCH(strspliti(names(DrugSensP),split=" ",i=1)))
abline(lm(CGTrtP~DrugSensP),lty=2)
legend('bottomleft',col='black',
       pch=getPCH(unique(strspliti(names(DrugSensP),split=" ",i=1))),
       legend=unique(strspliti(names(DrugSensP),split=" ",i=1)))
corVals <- cor.test(DrugSensP,CGTrtP)

title(sprintf('Change in pattern 2 vs cell survival\ncor %0.2f (p %0.2f)',
              corVals$estimate,corVals$p.value))

dev.off()


pdf('graphs/Pattern3_DrugSens.pdf',width=3.5,height=3.5)
CGTrtP <- CoGAPS_SDFIVE_3$Pmean[3,grep('Control',colnames(CoGAPS_SDFIVE_3$Pmean),invert=T)]

DrugSensP <- rep(NA, length(CGTrtP))
names(DrugSensP) <- names(CGTrtP)
for (i in names(DrugSensP)) {
  DrugSensP[i] <- DrugSens[strsplit(i,split=" ")[[1]][1],
                           strsplit(i,split=" ")[[1]][3]]
  CGTrtP[i] <- -CoGAPS_SDFIVE_3$Pmean[3,sub(strspliti(i,split=" ",i=3),
                                           'Control',i)] + 
    CGTrtP[i] 
}

plot(DrugSensP, CGTrtP, xlab='Cell Survival (%)',
     ylab='Pattern 3 in treatment - control',
     col=getCol(strspliti(names(DrugSensP),split=" ",i=3)),
     pch=getPCH(strspliti(names(DrugSensP),split=" ",i=1)))
abline(lm(CGTrtP~DrugSensP),lty=2)
legend('bottomleft',col='black',
       pch=getPCH(unique(strspliti(names(DrugSensP),split=" ",i=1))),
       legend=unique(strspliti(names(DrugSensP),split=" ",i=1)))
corVals <- cor.test(DrugSensP,CGTrtP)

title(sprintf('Change in pattern 3 vs cell survival\ncor %0.2f (p %0.2f)',
              corVals$estimate,corVals$p.value))

dev.off()

