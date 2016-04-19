sva.class2Model <- function(classes) {
  return(model.matrix(~factor(classes)))
}

### call to pSVA algorithm
psva <- function(dat, batch, ...) {
  
  # convert input class / categorical variables to a standard model matrix 
  if (class(batch)=='factor' | class(batch)=='character') {
    mod <- sva.class2Model(batch) 
  } else {
    stop('Invalid batch type for psva (require factor or character):', 
         class(batch))
  } 
  
  # find SV's
  psva.SV <- sva(dat=dat, mod=mod, ...)
  colnames(psva.SV$sv) <- paste('sv',1:ncol(psva.SV$sv))
  
  # fit data
  psva.fit <- lmFit(dat, cbind(mod,psva.SV$sv))
  
  # batch corrected data 
  psva.D <- sweep(psva.fit$coefficients[,paste('sv',1:ncol(psva.SV$sv))]%*%
                    t(psva.SV$sv), 1, psva.fit$coefficients[,"(Intercept)"], FUN="+")
  
  colnames(psva.D) <- colnames(dat)
  
  return(psva.D)
  
}