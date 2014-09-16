#' @title A-efficiency 
#' 
#' @description
#' \code{A_efficiency} A-efficiency of an existing design
#' @details
#' A-efficiencies for nested blocks design in a data frame with a
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
#' A_efficiency(blocks(treatments=50,replicates=4,blocklevels=c(4,5))$Design)
#' 
#' @export
A_efficiency=function(Design){
	strata=ncol(Design)-1
	nblocks=as.numeric(sapply(Design,nlevels))[1:strata]
	bounds=rep(0,strata)
	r=1/sqrt(tabulate(Design$Treatments))
	aeff=c(rep(0,strata))
	for (i in 1:strata) {
		k=1/sqrt(tabulate(Design[,i]))
		U=crossprod(t(crossprod(diag(r,nrow = length(r)),table(Design$Treatments,Design[,i]))),diag(k,nrow = length(k)))
		A=diag(length(r))-crossprod(t(U))
		aeff[i]=1/mean(1/eigen(A, symmetric=TRUE, only.values = TRUE)$values[1:length(r)-1])
		bounds[i]=upper_bounds(nrow(Design),nlevels(Design[,strata+1]),nblocks[i]) 
	}
	effic_dataframe=as.data.frame(cbind(nblocks, aeff, bounds))
	rnames=c("Main")
	if (strata>1)
	for (i in 1 : (strata-1)) rnames=c(rnames,paste("Sub",i))
	colnames(effic_dataframe)=c("Blocks","A-Efficiencies", "Upper Bounds")
	rownames(effic_dataframe)=rnames
	effic_dataframe	
} 
