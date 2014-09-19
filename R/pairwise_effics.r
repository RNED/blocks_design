#' @title Pairwise efficiencies
#' 
#' @description
#' \code{pairwise_efficiencies} Pairwise efficiencies of an existing design
#' @details
#'  Pairwise efficiencies for nested blocks design in a data frame with a
#'  column of block factor levels for the total blocks in each 
#'  nested blocks stratum and a final column containing the treatment 
#'  factor levels. The function returns a data frame with a row for 
#'  each stratum showing the total blocks, the A-efficiency factor and the 
#'  upper bound where available for each stratum of the design
#' 
#' @param Design A block design data frame such as returned by \code{\link{blocks}}
#' 
#' @examples 
#' 
#' # 4 replicates of 50 treatments in complete randomized blocks 
#' 
#' pairwise_efficiencies(blocks(treatments=50,replicates=4,blocklevels=c(4,5))$Design)
#' @export
pairwise_efficiencies=function(Design){
  nunits=nrow(Design)	
  strata=ncol(Design)-1
  pairwise_effics=vector(mode="list",length=strata)	
  TF=as.factor(Design$Treatments)
  T=matrix(0,nrow=nunits,ncol=nlevels(TF))	
  T[cbind(rep(1:nunits),TF)]=1
  r=apply(T, 2, sum)		
  for (i in 1:strata) {
    BF=as.factor(Design[,i])
    B=matrix(0,nrow=nunits,ncol=nlevels(BF))
    B[cbind(rep(1:nunits),BF)]=1
    B=B[ , (1:(ncol(B)-1)) ]
    B=crossprod(t(B),diag(1/sqrt(apply(B, 2, sum))))
    V=solve(diag(r,nrow = length(r))-crossprod(crossprod(B,T)))
    D= crossprod( t(rep(1,nlevels(TF)))  ,t(diag(V)))   +  crossprod( t(diag(V)), t(rep(1,nlevels(TF)))) - 2*V		
    N= crossprod( t(rep(1,nlevels(TF)))  ,1/t(r))       +  crossprod( 1/t(r), t(rep(1,nlevels(TF))))
    E=N/D
    E[upper.tri(E,diag=TRUE)] = NA
    pairwise_effics[[i]]=E
  }
  pairwise_effics
} 
