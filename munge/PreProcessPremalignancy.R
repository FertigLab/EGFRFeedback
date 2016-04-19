### fRMA normalize the Premalignancy data
Premalignancy.FRMA <- frma(Premalignancies.CEL)
Premalignancy.FRMA <- exprs(Premalignancy.FRMA)

# save frma data in the cache
ProjectTemplate::cache('Premalignancy.FRMA')

### sample annotation
sampleAnnot <- data.frame(coreID=c('FuDa'),
                          HAND.ID=c('FaDu'))
row.names(sampleAnnot) <- paste('JPer','Pre',sampleAnnot$coreID,'1a','U133Plus2_(HG-U133_Plus_2).CEL',sep='-')

Premalignancy.sampleAnnot <- sampleAnnot

# cache the data
ProjectTemplate::cache('Premalignancy.sampleAnnot')
