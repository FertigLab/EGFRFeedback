### Read in CEL files only for subsequent analyses to be performed
### during the munging steps (since cannot access data dir)

### libraries
library('affy')

### Premalignancies
Premalignancies.CEL <- ReadAffy(celfile.path='data/JPer050313/JPer050313 raw data/')
save(Premalignancies.CEL,file='data/Premalignancies.CEL.Rda')

### HaCaT
HaCaT.CEL.BATCH1 <- ReadAffy(celfile.path='data/JPer050313-2/JPer050313-2 raw data/')
save(HaCaT.CEL.BATCH1,file='data/HaCaT.CEL.BATCH1.Rda')

### HaCaT batch 2
HaCaT.CEL.BATCH2 <- ReadAffy(celfile.path='data/EFer080913/EFer080913 Raw Data/')
save(HaCaT.CEL.BATCH2,file='data/HaCaT.CEL.BATCH2.Rda')

### HaCaT batch 3
HaCaT.CEL.BATCH3 <- ReadAffy(celfile.path='data/MTha092214/MTha092214 raw data')
save(HaCaT.CEL.BATCH3,file='data/HaCaT.CEL.BATCH3.Rda')