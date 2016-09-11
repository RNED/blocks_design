#' @title Block designs
#'
#' @description
#' 
#' @details
#' 
#' @return
#' 
#' @references
#'
#' @examples
#'
#' @export
#' @importFrom stats anova lm
#'
testblocks = function(TF,replicates=1,order=NULL) {
  
  options(contrasts=c('contr.sum','contr.poly'))
  if (!is.data.frame(TF)) stop("Treatment factors not in a valid data frame") else print("Valid data frame")
  
  fInd=sapply(TF, is.factor)
  fTF=TF[, fInd,drop=FALSE]
  pTF=TF[,!fInd,drop=FALSE]
  
  if (ncol(fTF)>0) 
    for (i in 1: min(order,ncol(fTF))) 
      combF=combn(ncol(fTF),i)
    
 # TM=model.matrix(  as.formula(paste(" ~ ", paste(colnames(TF), collapse= "*")))  ,TF)
  
  colN=colnames(TF)
  o0=colN[c(1,2)]
  o1=colN[c(1,2)]
  form0=paste(o0, collapse= "+")
  form1=paste(o1, collapse= "*")
  TM=model.matrix(  as.formula( paste(" ~ ", form1 ) )  ,TF)
  #TM=model.matrix(  as.formula(paste(" ~ ", o1, collapse="+"))  ,TF)
  
  #TM=model.matrix(  as.formula(paste(" ~ ", paste(c(form0,form1), collapse="+")))  ,TF)
  

  if (ncol(pTF)>0) 
    for (i in 1: min(order,ncol(pTF))) {
      combP=combn(ncol(pTF),i)
      P=lapply(1:ncol(pTF), function(r){ poly(pTF[,r],degree=min(order,length(unique(pTF[,r]))-1),raw=TRUE)})
    }
  levelsF=sapply(fTF, nlevels)
  #model.matrix(~ a * poly(TFP,degree=4) , TF)
 
#TF=data.frame(a = factor(gl(2,25)), b = rep(rep(c(1,3,4,5,7),each=5),2),c = rep(c(1:5),10))
#TF=data.frame(b = rep(rep(c(1,3,4,5,7),each=5),2),c = rep(c(1:5),10))
#tempDesign=data.frame( do.call(cbind,lapply(1:ncol(pTF), function(r){ sample(nlevels(fDesign[,r]))[fDesign[,r]]})),seq_len(nrow(fDesign))) 

#Design=NULL
#for (i in 2:ncol(fTF)) {
#comb=combn(ncol(fTF),i)
#Design=cbind(Design,   do.call(cbind,lapply(1:ncol(comb), function(r){ apply(fTF[,comb[,r]],1,FUN="prod") }))  )
#}
#Design=cbind(fTF,Design)

#Design=scale( data.frame( rep(1:2,each=32), rep(rep(1:4,each=8),2), rep(rep(1:4,each=2),8), rep(1:2,32)) , center = TRUE, scale = FALSE)
#for (i in 2:4) {
 # comb=combn(4,i)
 # Design=cbind(Design,   do.call(cbind,lapply(1:ncol(comb), function(r){ apply(Design[,comb[,r]],1,FUN="prod") }))  )
#}
#colnames(Design)  <- 1:ncol(Design)

#comb=combn(ncol(fTF),2)

#do.call(cbind,lapply(1:ncol(comb), function(r){ apply(fTF[,comb[,r]],1,FUN="prod") }))

#Design=NULL
#for (i in 1:ncol(fTF)) {
#  comb=combn(ncol(fTF),i)
 # Design=cbind(Design,do.call(cbind,lapply(1:ncol(comb), function(r){  })))
#}
TM

}
#TF1=data.frame(gl(2,36), gl(3,12,72), gl(3,4,72),gl(4,1,72))
#TF2=data.frame(gl(2,36), gl(3,12,72), rep(rep(1:3,each=4),6), rep(1:4,18))
#testblocks(fTF)
#Q=qr(design)
#Q$rank