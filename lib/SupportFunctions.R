### function to get ith element from a string split
strspliti <- function (x, split, i = 1, ...) 
{
    return(geti(strsplit(x, split = split, ...), i))
}

### function to get the ith element from a list
geti <- function (l, i = 1) 
{
    if (class(l) != "list") {
        stop("geti cannot extract elements from objects of class %s", 
			 class(l))
    }
    iVal <- rep(NA, length(l))
    names(iVal) <- names(l)
    sidx <- which(sapply(l, length) >= i)
    if (length(sidx) > 0) {
        iVal[sidx] <- sapply(l[sidx], function(x) {
							 x[[i]]
							 })
    }
    return(iVal)
}

# get colors for each datapoint
getCol <- function(trt) {
  return(sapply(trt, 
                switch, 
                Control='black',
                Cetuximab='blue',
                Gefitinib='red',
                Afatinib='green'))
}

# get pch for each datapoint
getPCH <- function(cell) {
	return(sapply(cell,
				  switch,
				  "HaCaT-Mock"=16,
				  "HaCaT-EGFR"=15,
				  "HaCaT-PIK3CA"=17,
				  "HaCaT-HRAS"=18))
}
