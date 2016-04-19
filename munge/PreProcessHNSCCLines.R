# Using FRMA normalized HNSCC cell line data, obtained from another project

# extracting the expression data
HNSCC.CellLines.FRMA <- exprs(HNSCCCellLinesExprs)
ProjectTemplate::cache('HNSCC.CellLines.FRMA')

# extracting the sample annotations
HNSCC.CellLines.sampleAnnot <- pData(HNSCCCellLinesExprs)
ProjectTemplate::cache('HNSCC.CellLines.sampleAnnot')


