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
testblocks = function(TF,replicates=1,model=model) {
  options(contrasts=c('contr.sum','contr.poly'))
  TM=model.matrix(  (model)  ,TF)
TM
}
TF=data.frame(F1=gl(2,36), F2=gl(3,12,72), V1=rep(rep(1:3,each=4),6), V2=rep(1:4,18))
model=" ~ F1*F2 + V1*V2 + I(V1^2) + I(V2^2) + F1:V1 + F1:V2 + F2:V1 + F2:V2"
TM=model.matrix(  as.formula(model)  ,TF)


#TM=model.matrix(  as.formula( paste(" ~ ", model  ) )  ,TF)

#m=model.matrix(~ F1*F2 + V1*V2 + I(V1^2) + I(V2^2) + F1:V1 + F1:V2 + F2:V1 + F2:V2,TF)
