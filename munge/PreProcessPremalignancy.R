### fRMA normalize the Premalignancy data
Premalignancy.FRMA <- frma(Premalignancies.CEL)
Premalignancy.FRMA <- exprs(Premalignancy.FRMA)

# save frma data in the cache
ProjectTemplate::cache('Premalignancy.FRMA')

### sample annotation
sampleAnnot <- data.frame(coreID=c(34692,37428,37435,37445,37448,37461,37480,'FuDa'),
                          HAND.ID=c('0000034672','0000037428','0000037435','0000037445','0000037448','0000037461','0000037480','FaDu'),
                          Patient=c('F','H','H','H','F','F','F',NA),
                          Treatment=c('Post','Post',rep('Pre',5),NA),
                          RIN=c(9.3,6.9,8.3,6.3,7.1,7.5,8,9.8))
row.names(sampleAnnot) <- paste('JPer','Pre',sampleAnnot$coreID,'1a','U133Plus2_(HG-U133_Plus_2).CEL',sep='-')

Premalignancy.sampleAnnot <- sampleAnnot

# cache the data
ProjectTemplate::cache('Premalignancy.sampleAnnot')
