#' @title Efficiencies of pairwise differences
#' 
#' @description
#' \code{efficiencies} finds efficiencies of pairwise treatment differences
#' 
#' @details
#'  Efficiencies of pairwise treatment differences for a design built by using the \code{blocks} function. 
#'  The function calculates a two-way table showing the efficiency of each pairwise treatment difference. 
#'  For equi-replicate designs only, the harmonic mean of the pairwise efficiency factors 
#'  should equal the A-efficiency factor calculated by the \code{blocks} function. 
#'  
#' @param Design A block design data frame as returned by \code{\link{blocks}}
#' 
#' @return  
#' \item{Incidences[[i]]}{Blocks by treatments incidence matrices  i=1... for each stratum in the design}
#' 
#' @examples 
#' 
#' # 4 replicates of 50 treatments in complete randomized blocks 
#' 
#' efficiencies(blocks(treatments=50,replicates=4,blocklevels=c(4,5))$Design)
#' @export
efficiencies=function(Design){
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
    B=B[ , (1:(ncol(B)-1)),drop=FALSE ]
    B=crossprod(t(B),diag(1/sqrt(apply(B, 2, sum)),nrow=ncol(B)))
    V=solve(diag(r,nrow = length(r))-crossprod(crossprod(B,T)))
    D= crossprod( t(rep(1,nlevels(TF)))  ,t(diag(V)))   +  crossprod( t(diag(V)), t(rep(1,nlevels(TF))) ) - 2*V		
    N= crossprod( t(rep(1,nlevels(TF)))  ,1/t(r))       +  crossprod( 1/t(r), t(rep(1,nlevels(TF))))
    E=N/D
    E[upper.tri(E,diag=TRUE)] = NA
    pairwise_effics[[i]]=E
  }
  pairwise_effics
} 
