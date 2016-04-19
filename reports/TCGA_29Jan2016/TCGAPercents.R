library('xlsx')

TCGA.file <- '~/Dropbox/K25_Fertig/Aim3Analysis/reports/TCGA_29Jan2016/cBioPortal_Alterations_1Dec2015.xlsx'

TCGA <- loadWorkbook(TCGA.file)
TCGA.sheets <- getSheets(TCGA)

TCGA.data <- lapply(names(TCGA.sheets),
                    function(x){
                      x <- read.xlsx(file=TCGA.file,sheetName=x,stringsAsFactors=F)
                      row.names(x) <- x$Case.ID
                      x <- x[,setdiff(sort(colnames(x)),'Case.ID')]
                      return(x)
                    })
names(TCGA.data) <- names(TCGA.sheets)

TCGA.data$HNSC <- rbind(TCGA.data$HNSC_HPVNeg,
                        TCGA.data$HNSC_HPVPos[,colnames(TCGA.data$HNSC_HPVNeg)])
TCGA.data <- TCGA.data[setdiff(names(TCGA.data),c('HNSC_HPVPos','HNSC_HPVNeg'))]

TCGA.data.all <- TCGA.data[[1]]
row.names(TCGA.data.all) <- paste(names(TCGA.data)[1],row.names(TCGA.data[[1]]),sep="-")
for (i in 2:length(TCGA.data)) {
  tmp <- TCGA.data[[i]]
  row.names(tmp) <- paste(names(TCGA.data)[i],
                          row.names(TCGA.data[[i]]),sep="-")
  
  TCGA.data.all <- rbind(TCGA.data.all,tmp[,colnames(TCGA.data.all)])
}

calcPercents <- function(TCGA.table) {
  genes <- list(EGFR='EGFR',
                PI3K=paste0('PIK3C',c('A','B','G','D')),
                RAS=paste0(c('H','K','N'),'RAS'),
                AKT=paste0('AKT',1:3),
                RAF=c('BRAF','ARAF'),
                STAT3='STAT3',
                STAT5B='STAT5B',
                NFKB=c('NFKB1','NFKB2','RELA','RELB'),
                FOXO=c('FOXO1'),
                MEK=c('MAP2K7'),
                ERK='MAPK1',
                ELK1='ELK1',
                MYC='MYC')
  return(sapply(genes,function(x){
    if (length(x)==1) {
      return(sum(!is.na(TCGA.table[,x]))/nrow(TCGA.table))
    } else {
      return(sum(apply(!is.na(TCGA.table[,x]),1,any))/nrow(TCGA.table))
    }
  }))
}

TCGA.percents <- sapply(TCGA.data,calcPercents)

write.table(cbind(row.names(TCGA.percents),TCGA.percents),
            file='~/Dropbox/K25_Fertig/Aim3Analysis/reports/TCGA_29Jan2016/TCGAPercents.csv',
            sep=",",row.names=F)

round(100*apply(TCGA.percents,1,mean))
