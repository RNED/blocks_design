#' @title Block designs 
#' 
#' @description
#' 
#' Constructs randomized nested block designs for unstructured treatment sets where treatments can have any arbitrary levels of replication
#' and blocks can have any arbitrary feasible depth of nesting.
#' 
#' @details
#' 
#' The \code{treatments} and \code{replicates} parameters partition the total number of treatments into sets of equally replicated treatments where 
#' \code{treatments} contains the set sizes and \code{replicates} contains the set replication numbers. 
#' The sum of the set sizes is the total number of treatments and the sum of the cross-products of the set sizes and the replication numbers
#' is the total number of plots. Treatments are numbered consecutively according to the ordering of the consecutive treatment sets. 
#' 
#' The \code{nestedLevels} parameter contains the nested blocks for each individual stratum taken in order from the highest to the lowest.
#' The first number is the number of main blocks, the second, if any, is the number of sub-blocks nested in each main block, the third, if any, 
#' is the number of sub-sub-blocks nested in each sub-block,and so on for all the reqired strata. If left blank, the default block design is a 
#'  maximal set of orthogonal main blocks, where the maximal number of of orthogonal main blocks is the highest common factor of the replication numbers. 
#'
#' The block sizes in each blocks stratum are always as equal as possible and never differ by more than a single unit. If the number of nested blocks 
#' in a particular stratum exactly divdes the number of units, the block sizes in that stratum will be exactly equal; otherwise the block sizes
#' will differ by at most one unit.
#' 
#'  Special square and rectangular lattice designs are constructed algebraically and include all designs that can be constructed from 
#'  a single latin square, from a complete sets of prime or prime-power orthogonal latin squares or from a pair of orthogonal 10 x 10 Latin squares. 
#'  All other non-orthogonal block designs are constructed by a D-optimality swapping algorithm that makes improving swaps between 
#'  blocks within the same stratum until a local optima is atttained. The swapping algorithm works from the top stratum downwards and
#'  is always constrained to make improving swaps within the levels of any existing blocks. The whole process is repeated according to the 
#'  number of searches defined by the search parameter and the design returned will be the design with the best overall stratum efficiencies in each stratum 
#'  taken in top-down order.
#'  
#'  Lattice designs where v is a prime-power require the \code{\link[crossdes]{MOLS}} package.
#' 
#'  The principle design outputs comprise:
#' \itemize{
#'  \item  A design matrix showing the allocation of treatments to blocks with successive nested blocks factors arranged in successive columns in standard block order.  \cr
#'  \item  A design matrix as above but with the last (bottom) blocks factor shown arranged horizontally to give a plan view. \cr
#'  \item  A set of incidence matrices, one for each blocks stratum, showing the number of times each treatment occurs in each block for each stratum. \cr
#'  \item  A table showing the achieved D- and A-efficiency factors for each nested blocks stratum together with an A-efficiency upper bound, where available. \cr
#'  \item  A table showing a skeleton analysis of degrees of freedom for the combined block and treatment design. \cr
#' } 
#' 
#' Very occasionally, the algorithm may fail to converge due to a near-singular design with a large number of single plot blocks.
#' In that case, it may be best to build a simpler design with larger blocks and then to add the extra block constraints by hand using ad hoc or heuristic methods.     
#' 
#' @param treatments treatment numbers giving a partition of the total required number of treatments into sets of equally replicated treatments.
#' 
#' @param replicates replication numbers giving the replictaion for each set of equally replicated treatments defined by the \code{treatments} partition.
#' 
#' @param nestedLevels factor levels that define the number of nested blocks in each succesive blocks stratum taken in order from the highest to the lowest. 
#' The default is the hcf of the replication numbers.
#' 
#' @param seed integer initializing the random number generator. The default is a random seed.
#' 
#' @param searches maximum number of local optima searched for a design optimization. The default is 1 plus the floor of 2000 divided by the number of model parameters.
#' 
#' @param jumps number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#' 
#' @return  
#' \item{Design}{Data frame giving the optimized block and treatment factors in plot order}
#' \item{Plan}{Data frame giving a plan view of the treatments in the bottom stratum of the design}
#' \item{AOV}{Data frame giving a skeleton analysis of variance of the degrees of freedom of the design}
#' \item{Incidences}{Blocks-by-treatments incidence matrices for each stratum of the design}
#' \item{Efficiencies}{The achieved A- and D-efficiencies for each stratum of the design together with an A-efficiency upper-bound, where available}
#' \item{seed}{Numerical seed for random number generator}
#' \item{searches}{Maximum number of searches in each stratum}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima}
#' 
#' @references
#' 
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. http://CRAN.R-project.org/package=crossdes
#' 
#' @examples
#' 
#' # 3 treatments x 2 replicates, 2 treatments x 4 replicates and 4 treatments x 3 replicates  
#' # the hcf of the replication numbers is 1 therefore the default design is completely randomized 
#' blocks(treatments=c(3,2,4),replicates=c(2,4,3))
#' 
#' # 4 treatments x 4 replicates with 2 main blocks each containing two complete replicates  
#' blocks(treatments=4,replicates=4,blocklevel=2)
#' 
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block 
#' blocks(treatments=50,replicates=4,nestedLevels=c(4,5))
#' 
#' # as above but with 20 additional single replicate treatments 
#' # giving exactly one single replicate treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),nestedLevels=c(4,5))
#' 
#' # 64 treatments x 2 replicates with 2 main blocks and five succesively nested 2-level factors
#' blocks(treatments=64,replicates=2,nestedLevels=c(2,2,2,2,2,2))
#' 
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,nestedLevels=4)
#' 
#' # concurrence matrix of balanced incomplete block design 
#' crossprod(blocks(13,4,13,searches=100)$Incidences[[1]])
#' 
#' # concurrence matrix for 13 treatments x 4 replicates and 13 treatments with one rep in 13 blocks 
#' crossprod(blocks(c(13,13),c(4,1),13)$Incidences[[1]])
#' 
#' # 2**10 treatments x 2 replicates in 2**10 blocks giving a fully saturated blocks design 
#' # (requires a considerable time to run!)
#' \dontrun{ d=blocks(1024,2,rep(2,10)) }
#'          
#' @export
#' @importFrom stats anova lm
#' 
blocks = function( treatments, replicates, rowLevels=HCF(replicates),colLevels=NULL,searches=(1+2000%/%(sum(treatments*replicates))),seed=sample(10000,1),jumps=1) { 
 
# ******************************************************************************************************************************************************** 
#  Generates a vector of block sizes for a particular stratum where all blocks are as equal as possible and never differ by more than a single unit
# ********************************************************************************************************************************************************
  Sizes=function(sizes,nestedLevels) { 
    for  (i in 1:length(nestedLevels)) {    
      newsizes=NULL
      for (z in 1: length(sizes)) 
        newsizes=c(newsizes, rep(sizes[z] %/% nestedLevels[i], nestedLevels[i]) + 
                     c( rep(1, sizes[z] %% nestedLevels[i]), rep(0,(nestedLevels[i]-sizes[z] %% nestedLevels[i])))) 
      sizes=newsizes
    }   
    sizes 
  } 
  
  # ******************************************************************************************************************************************************** 
  #  Generates a vector of block sizes for a latin square
  # ********************************************************************************************************************************************************
  BlockSizes=function(blocklevs,nunits) {
    rblocks=nunits%/%blocklevs[1]+c(rep(1,nunits%%blocklevs[1]),rep(0,(blocklevs[1]-nunits%%blocklevs[1])))
    sblocks = vector(length=blocklevs[1]*blocklevs[2])	
    shift = 0
    for (j in 1: blocklevs[1]) {
      bresid = rblocks[j]%%blocklevs[2];
      bbase = rblocks[j]%/%blocklevs[2];
      for (k in 1: blocklevs[2])
        sblocks[ (j-1)*blocklevs[2] + k]= bbase
      if (bresid > 0) {
        for (z in 0: (bresid-1))
          sblocks[(j-1)*blocklevs[2] + (z + shift)%%blocklevs[2] +1]=bbase+1
        shift = shift+bresid;
      }
    }
    sblocks
  }
  
# ******************************************************************************************************************************************************** 
# Finds the highest common factor (hcf) of a set of numbers omitting any zero values (Euclidean algorithm)
# ********************************************************************************************************************************************************
  HCF=function(replevs)  {
    replevs=sort(replevs[replevs>0])
    for (i in 1: length(replevs))
      while (!isTRUE(all.equal(replevs[i]%%replevs[1],0))) replevs[c(1,i)] = c(replevs[i]%%replevs[1], replevs[1])
      replevs[1]
  }   
# ******************************************************************************************************************************************************** 
# Tests a given number for primality and returns TRUE or FALSE
# ********************************************************************************************************************************************************
  isPrime=function(v) {
    if (v < 4) return(TRUE) 
    if ( isTRUE(all.equal(v %% 2,0)) ||  isTRUE(all.equal(v %% 3,0)) ) return(FALSE) 
    if (v<25) return(TRUE)
    for(i in  6*rep(1:floor((sqrt(v)+1)/6)) )
      if ( isTRUE(all.equal(v %% (i-1) , 0)) ||   isTRUE(all.equal(v %% (i+1) , 0)) ) return(FALSE) 
    return(TRUE)
  } 
# ******************************************************************************************************************************************************** 
# Contrasts for factor NF centered within the levels of factor MF to ensure that NF information is estimated within the levels of factor MF only  
# ********************************************************************************************************************************************************
 Contrasts=function(MF,NF) {
    NM=matrix(0,nrow=length(NF),ncol=nlevels(NF))
    NM[cbind(1:length(NF),NF)]=1 # factor indicator matrix  
    for (i in 1:nlevels(MF)) 
      NM[MF==i,]=scale(NM[MF==i,] , center = TRUE, scale = FALSE)
  NM
  }
# ******************************************************************************************************************************************************** 
# Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
# mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
# 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb   
# ********************************************************************************************************************************************************
  UpDate=function(MTT,MBB,MTB,ti,tj,bi,bj) {  
    mtt=MTT[ti,ti]+MTT[tj,tj]-2*MTT[tj,ti]
    mbb=MBB[bi,bi]+MBB[bj,bj]-2*MBB[bi,bj]
    mtb=1-MTB[ti,bi]+MTB[tj,bi]+MTB[ti,bj]-MTB[tj,bj]  
    TBbij=MTB[,bi]-MTB[,bj]
    TBtij=MTB[ti,]-MTB[tj,]
    TTtij=MTT[,ti]-MTT[,tj]
    BBbij=MBB[bi,]-MBB[bj,]
    Z1 = (TBbij-TTtij)/sqrt(2*mtb+mtt+mbb)   
    Z2 = (BBbij-TBtij)/sqrt(2*mtb+mtt+mbb)
    W1 = ( sqrt(2*mtb+mtt+mbb)*(TTtij+TBbij) - (mbb-mtt)*Z1) /(2*sqrt(mtb**2-mtt*mbb))
    W2 = ( sqrt(2*mtb+mtt+mbb)*(TBtij+BBbij) - (mbb-mtt)*Z2) /(2*sqrt(mtb**2-mtt*mbb))
    MTT = MTT - tcrossprod(Z1) + tcrossprod(W1)
    MBB = MBB - tcrossprod(Z2) + tcrossprod(W2)
    MTB = MTB - tcrossprod(Z1,Z2) + tcrossprod(W1,W2) 
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }  
  # ******************************************************************************************************************************************************** 
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  optEffics=function(TF,BF) { 
    r=nlevels(TF)
    k=nlevels(BF)
    if (r<=k) 
      e=eigen( (diag(r)-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(r-1)] else    
      e=c(rep(1,(r-k)),
          eigen((diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])  
    round(c(mean(e)*prod(e/mean(e))^(1/length(e)),1/mean(1/e)),6)
  }
# ******************************************************************************************************************************************************** 
# Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
# Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
# ********************************************************************************************************************************************************
DMax=function(MTT,MBB,MTB,TF,BF,Restrict){   
  relD=1
  mainSizes=tabulate(Restrict)
  nSamp=pmin(rep(8,nlevels(Restrict)),mainSizes)
  repeat {
    improved=FALSE
    for (k in 1:nlevels(Restrict)) {
      S=sort(sample((1:length(TF))[Restrict==k],nSamp[k])) 
      TB=MTB[TF[S],BF[S],drop=FALSE]-tcrossprod(MTB[cbind(TF[S],BF[S])],rep(1,nSamp[k]))
      dMat=(TB+t(TB)+1)**2-
        (2*MTT[TF[S],TF[S],drop=FALSE]-tcrossprod(MTT[cbind(TF[S],TF[S])]+rep(1,nSamp[k]) ) + tcrossprod(MTT[cbind(TF[S],TF[S])]) + 1)*
            (2*MBB[BF[S],BF[S],drop=FALSE]-tcrossprod(MBB[cbind(BF[S],BF[S])]+rep(1,nSamp[k]) ) + tcrossprod(MBB[cbind(BF[S],BF[S])]) + 1)
      sampn=which.max(dMat) 
      i=1+(sampn-1)%%nSamp[k]
      j=1+(sampn-1)%/%nSamp[k]
      if ( !isTRUE(all.equal(dMat[i,j],1)) && dMat[i,j]>1) {
        improved=TRUE
        relD=relD*dMat[i,j]
        up=UpDate(MTT,MBB,MTB,TF[S[i]],TF[S[j]],BF[S[i]],BF[S[j]])
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        TF[c(S[i],S[j])]=TF[c(S[j],S[i])]
      }
    } 
    if (improved) next
    if (sum(nSamp) < min(length(TF),512)) nSamp=pmin(mainSizes,2*nSamp) else break
  }  
  list(MTT=MTT,MBB=MBB,MTB=MTB,TF=TF,relD=relD)
}  
# ******************************************************************************************************************************************************** 
#  Number of searches for an optimization with selected number of searches and selected number of junps to escape local optima
# ********************************************************************************************************************************************************
  Optimise=function(TF,BF,MTT,MBB,MTB,searches,jumps,Restrict)  {
    globrelD=0
    relD=1
    globTF=TF
    treps=tabulate(TF)
    breps=tabulate(BF)
    if (identical(max(treps),min(treps)) && identical(max(breps),min(breps))  )
      bound=upper_bounds(length(TF),nlevels(TF),nlevels(BF)) else bound=NA
    for (r in 1 : searches) {
      dmax=DMax(MTT,MBB,MTB,TF,BF,Restrict)
      if ( !isTRUE(all.equal(dmax$relD,1)) && dmax$relD>1) {
        relD=relD*dmax$relD
        TF=dmax$TF
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB 
        if (!isTRUE(all.equal(relD,globrelD)) && relD>globrelD) {
          globTF=TF
          globrelD=relD
          if ( !is.na(bound) && isTRUE(all.equal(bound,optEffics(globTF,BF)[2]))) break
        }
      }
      if (r==searches) break
      for (iswap in 1 : jumps) {
        counter=0
        repeat { 
          counter=counter+1
          s1=sample(1:length(TF),1)
          z=(1:length(TF))[Restrict==Restrict[s1] & BF!=BF[s1] & TF!=TF[s1]]
          if (length(z)==0) next
          if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z[1])
          dswap = (1+MTB[TF[s[1]],BF[s[2]]]+MTB[TF[s[2]],BF[s[1]]]-MTB[TF[s[1]],BF[s[1]]]-MTB[TF[s[2]],BF[s[2]]])**2-
            (2*MTT[TF[s[1]],TF[s[2]]]-MTT[TF[s[1]],TF[s[1]]]-MTT[TF[s[2]],TF[s[2]]])*(2*MBB[BF[s[1]],BF[s[2]]]-MBB[BF[s[1]],BF[s[1]]]-MBB[BF[s[2]],BF[s[2]]])  
          if (dswap>.1 | counter>1000) break
        }
        if (counter>500) return(globTF) # no non-singular swaps
        relD=relD*dswap
        up=UpDate(MTT,MBB,MTB,TF[s[1]],TF[s[2]], BF[s[1]], BF[s[2]])
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        TF[c(s[1],s[2])]=TF[c(s[2],s[1])]  
      } 
    }
   globTF
  } 
  # ******************************************************************************************************************************************************** 
  # Random swaps
  # ********************************************************************************************************************************************************    
  Swaps=function(TF,MF,BF,pivot,rank) {
    n=length(TF)
    candidates=NULL
    while (isTRUE(all.equal(length(candidates),0))) {
      if (rank<(n-1)) s1=sample(pivot[(1+rank):n],1) else s1=pivot[n]
      candidates = rep(1:n)[MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1]]
    }
    if ( length(candidates)>1 )
        s2=sample(candidates,1) else s2=candidates[1] 
    s=c(s1,s2)
  }
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  NonSingular=function(TF,MF,BF) { 
      BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
      BM[cbind(1:length(BF),BF)]=1
      fullrank=nlevels(TF)+nlevels(BF)-1
      TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
      TM[cbind(1:length(TF),TF)]=1
      Q=qr(t(cbind(BM,TM)))
      rank=Q$rank
      pivot=Q$pivot
      times=0
      while (rank<fullrank & times<1000) {
        times=times+1
        s=Swaps(TF,MF,BF,pivot,rank)
        rindex=(1:length(TF))
        rindex[c(s[1],s[2])]=rindex[c(s[2],s[1])]
        newQ=qr(t(cbind(BM,TM[rindex,])))
        if (isTRUE(all.equal(newQ$rank,rank)) || newQ$rank>rank) { 
          TF=TF[rindex]
          TM=TM[rindex,]
          rank=newQ$rank
          pivot=newQ$pivot
        } 
      }
      if ( isTRUE( times > 999 ) )
         stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
    return(TF)
  }  
  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  GenOpt=function(TF,Design,searches,jumps,stratum,nestedLevels,hcf) { 
    Design=cbind(as.factor(rep(1,nrow(Design))),Design)
    MF=Design[,stratum]
    BF=Design[,stratum+1]
    rand=sample(1:length(TF))
    TF=TF[rand][order(MF[rand])]
    if (!isTRUE(all.equal(hcf %% prod(nestedLevels[1:stratum]),0))) 
    TF=NonSingular(TF,MF,BF)
    if (length(TF)==0) return(TF)
    blevels=nlevels(BF)%/%nlevels(MF)
    BM=Contrasts(MF,BF)[, rep(c(rep(TRUE,(blevels-1)),FALSE),nlevels(MF)),drop=FALSE]
    TM=Contrasts(MF,TF)[,-nlevels(TF),drop=FALSE] 
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    MBB=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
    MTT=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))  
    MTB=matrix(0,nrow=nlevels(TF),ncol=nlevels(BF))
    MBB[1:(nlevels(BF)-nlevels(MF)), 1:(nlevels(BF)-nlevels(MF))]=V[1:(nlevels(BF)-nlevels(MF)),1:(nlevels(BF)-nlevels(MF)),drop=FALSE]
    MTT[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[(nlevels(BF)-nlevels(MF)+1):ncol(V),(nlevels(BF)-nlevels(MF)+1):ncol(V), drop=FALSE]
    MTB[1:(nlevels(TF)-1),  1:(nlevels(BF)-nlevels(MF)) ]=V[(nlevels(BF)-nlevels(MF)+1):ncol(V),1:(nlevels(BF)-nlevels(MF)),drop=FALSE]
    perm=order(order( (1:nlevels(BF))%%blevels ==0  ))  
    MTB=MTB[,perm]
    MBB=MBB[perm,perm] 
    TF=Optimise(TF,BF,MTT,MBB,MTB,searches,jumps,MF)
    TF
  }  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  RC_Opt=function(Design,searches,jumps,strata) { 
    strata=(ncol(Design)-1)/2
    TF=as.factor(Design[,ncol(Design)])
    NestInt =data.frame(rep(1,nrow(Design)), do.call(cbind,lapply(1:strata, function(r){ (as.numeric(Design[,2*r-1])-1)*nlevels(Design[,2*r])+as.numeric(Design[,2*r])})))
    for (i in 1:strata)
    NestInt[,i+1] = (NestInt[,i]-1)*max(NestInt[,i+1])+NestInt[,i+1]
    NestInt[]=lapply(NestInt, as.factor)
    NestRows=data.frame(rep(1,nrow(Design)))
    for (i in 1:strata)
      NestRows=cbind( NestRows, (as.numeric(NestInt[,i])-1)*nlevels(Design[,2*i-1])+as.numeric(Design[,2*i-1]))
    NestRows[]=lapply(NestRows, as.factor)
    NestCols=data.frame(rep(1,nrow(Design)))
    for (i in 1:strata)
      NestCols=cbind( NestCols, (as.numeric(NestInt[,i])-1)*nlevels(Design[,2*i])+as.numeric(Design[,2*i]))
    NestCols[]=lapply(NestCols, as.factor) 
    for (i in 1:strata){
      MF=NestInt[,i]
        for (j in 1:2) {
        if (j==1) {
          BF=NestRows[,i+1]
          Res=MF
        } else {
          BF=NestCols[,i+1]
          Res=NestRows[,i+1]
        }
        blevels=nlevels(BF)%/%nlevels(MF)
        if (blevels<2) next
        TF=NonSingular(TF,Res,BF)
        TM=Contrasts(MF,TF)[,-nlevels(TF),drop=FALSE] 
        BM=Contrasts(MF,BF)[, rep(c(rep(TRUE,(blevels-1)),FALSE),nlevels(MF)),drop=FALSE]
        V=chol2inv(chol(crossprod(cbind(BM,TM))))
        MBB=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
        MTT=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))  
        MTB=matrix(0,nrow=nlevels(TF),ncol=nlevels(BF))
        MBB[1:(nlevels(BF)-nlevels(MF)), 1:(nlevels(BF)-nlevels(MF))]=V[1:(nlevels(BF)-nlevels(MF)),1:(nlevels(BF)-nlevels(MF)),drop=FALSE]
        MTT[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[(nlevels(BF)-nlevels(MF)+1):ncol(V),(nlevels(BF)-nlevels(MF)+1):ncol(V), drop=FALSE]
        MTB[1:(nlevels(TF)-1),1:(nlevels(BF)-nlevels(MF))]=V[(nlevels(BF)-nlevels(MF)+1):ncol(V),1:(nlevels(BF)-nlevels(MF)),drop=FALSE]
        perm=order(order((1:nlevels(BF))%%blevels ==0))  
        MTB=MTB[,perm]
        MBB=MBB[perm,perm] 
        TF=Optimise(TF,BF,MTT,MBB,MTB,searches,jumps,Res)
        }
    }
  Design[,ncol(Design)]=as.factor(TF)
  return(Design)
  }  
# ******************************************************************************************************************************************************** 
# Generates an initial orthogonal design then builds algebraic lattice blocks or calls the general block design algorithm as appropriate
# ********************************************************************************************************************************************************     
    optTF=function(Design,treatlevs,replevs,nestedLevels,searches,jumps) {
    nunits=sum(treatlevs*replevs)
    ntrts=sum(treatlevs)
    hcf=HCF(replevs)
    MB=rep(1:hcf,each=(nunits%/%hcf))
      repeat{
      TF=rep( rep(1:ntrts,rep(replevs%/%hcf,treatlevs)), hcf)
      rand=sample(nunits)
      TF=as.factor( TF[rand][order(MB[rand])] )
      firstPass=TRUE
      v=sqrt(ntrts)  
      regrep=isTRUE(all.equal(max(replevs),min(replevs)))
      for (stratum in 1 : length(nestedLevels)) { 
        nblocks=prod(nestedLevels[1:stratum])
        if (isTRUE(all.equal(hcf%%nblocks,0))) next
        regular=firstPass &&  regrep && isTRUE(all.equal(nunits%%nblocks,0))
        firstPass=FALSE
        sqrLattice  = regular && isTRUE(all.equal(v,floor(v))) && isTRUE(all.equal(nunits,v*nblocks))
        w=sqrt(nblocks)
        s=ntrts/w
        rectLattice = regular && isTRUE(all.equal(replevs[1],w)) && isTRUE(all.equal(s,floor(s))) && (s<w)
        if (sqrLattice  && replevs[1]<4) {
            t=c(rep(0:(v-1),each=v),rep(0:(v-1),v)+v)
            if (replevs[1]>2)
            for (i in 0: (v-1)) 
              t=c(t,(rep(0:(v-1))+i)%%v + 2*v)
            TF=as.factor(rep(1:(v*v),replevs[1])[order(t)])
        } else if (sqrLattice  &&  replevs[1]<(v+2)  && isPrime(v) ) {
            t=c(rep(0:(v-1),each=v),rep(0:(v-1),v)+v)
            for (z in 1: (replevs[1]-2))
              for (j in 0: (v-1)) 
                t=c(t,(rep(0:(v-1))*z +j)%%v + v*(z+1) )
              TF=as.factor(rep(1:(v*v),replevs[1])[order(t)])
        } else if (sqrLattice  && replevs[1]<(v+2)  &&  ntrts%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
              index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==ntrts)
              mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])			
              TF=c(rep(1:ntrts), rep(1:ntrts)[order(rep(0:(v-1),v))])
              for (i in 1: (replevs[1]-2))
                TF=c(TF, rep(1:ntrts)[order(as.numeric(mols[,,i]) ) ]) 
              TF=as.factor(TF)
        } else if (sqrLattice  && v==10  && replevs[1]==4) {
          TF=as.factor(c(
            rep(1:100), rep(1:100)[order(rep(0:9,10))],
            rep(1:100)[order(c(
              1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
              5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8))],
            rep(1:100)[order(c(
              1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
              6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6))]
          )) 
        } else if (rectLattice  && isPrime(w) ) {
              for (z in 1:(nunits/nblocks))
                for (j in 0: (w-1)) 
                  for (k in 0: (w-1)) 
                    TF[ (k + j*w)*nunits/nblocks + z]=(j+k*z)%%w + (z-1)*w +1
                  TF=as.factor(TF)
        } else if (rectLattice &&  nblocks%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
              index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==nblocks)
              mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])		
              TF=NULL
              for (i in 1: w)
                for (j in 1: w)
                  for (k in 1: (nunits/nblocks))
                    TF=c(TF,mols[i,j,k]+(k-1)*w)
              TF=as.factor(TF)
        } else TF=GenOpt(TF,Design,searches,jumps,stratum,nestedLevels,hcf)
      }
      if (length(TF)>0) break
    }
  levels(TF)=sample(1:ntrts) 
  TF
  }
# ******************************************************************************************************************************************************** 
# Finds efficiency factors for the blocks in each stratum of a design 
# ********************************************************************************************************************************************************     
  A_Efficiencies=function(Design,treatments,replicates,rowLevels,colLevels)     {
    strata=length(rowLevels)
    nestLevels=rowLevels*colLevels
    cumnestLevs=c(1,cumprod(nestLevels))
    hcf=HCF(replicates)
    nunits=nrow(Design)
    bounds=rep(NA,2*strata) 
    effics=matrix(nrow=2*strata,ncol=2)
    
    for (i in 1:strata) { 
      if (isTRUE(all.equal(max(replicates),min(replicates))) ) {
        nestrows=cumnestLevs[i]*rowLevels[i]
        if (nunits%%nestrows==0) {
          bounds[2*i-1]=upper_bounds(nunits,sum(treatments), nestrows) 
          } else {
          bounds[i]=1
          }
        
        nestcols=cumnestLevs[i]*colLevels[i]
        if (nunits%%nestcols==0) {
          bounds[2*i]=upper_bounds(nunits,sum(treatments), nestcols) 
        } else {
          bounds[i]=1
        }
      }
        effics[2*i-1,]=optEffics(Design$Treatments,Design[,2*i-1])
        effics[2*i,]=optEffics(Design$Treatments,Design[,2*i])
    }

    blocks=NULL
    for (i in 1:strata)
    blocks=c(blocks, rbind(    nlevels(Design[,2*i-1]) , nlevels(Design[,2*i])            )              )

    stratanames=rep("Main",2)
    if (strata>1)
    for (i in 1:(strata-1)) stratanames=c(stratanames,rep(paste0("Sub_",i),2))
    
    rowscols=rep(c("Rows","Cols"),strata)
    
    efficiencies=data.frame(cbind(stratanames,rowscols,blocks, effics, bounds)) 
    
    colnames(efficiencies)=c("Stratum","Rows:Cols","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies[,'Blocks'] = as.factor(efficiencies[,'Blocks'])
    efficiencies
  }
# ******************************************************************************************************************************************************** 
# Carries out some input validation
# ********************************************************************************************************************************************************     
 testInputs=function(treatments,replicates,rowLevels,colLevels,searches,seed,jumps) {  
   if (missing(treatments) | missing(replicates) )  
     return(" Treatments or replicates not defined ")   
   if (is.null(treatments) | is.null(replicates))  
     return(" Treatments or replicates list is empty ")   
   if (anyNA(treatments) | anyNA(replicates) ) 
     return(" NA values not allowed")
   if (any(!is.finite(treatments)) | any(!is.finite(replicates)) | any(is.nan(treatments)) | any(is.nan(replicates))) 
     return(" Treatments and replicates can contain only finite integers ")
   if ( length(treatments)!=length(replicates) ) 
     return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))
  if (any(treatments<1)) 
     return("Treatments must be non-negative integers")
  if (any(replicates<1)) 
    return("Replicates must be non-negative integers")  
   if (!is.null(rowLevels)) {
     if (anyNA(rowLevels) ) return(" NA rowLevels values not allowed") 
     if (any(!is.finite(rowLevels)) | any(is.nan(rowLevels)) ) return(" rowLevels can contain only finite integers ")
     if (min(rowLevels)<1) return (" rowLevels must be at least one ")
   }
   if (!is.null(searches)) {
     if (anyNA(searches) ) return(" NA searches values not allowed") 
     if ( any(!is.finite(searches)) | any(is.nan(searches))) return(" Searches must be a finite integer ") 
     if (searches<1)  return(" Repeats must be at least one ")   
   }  
  if (!is.null(jumps)) {
    if (anyNA(jumps) ) return(" NA jumps values not allowed") 
    if ( !all(is.finite(jumps)) | !all(!is.nan(jumps))) return(" jumps must be a finite integer ") 
    if (jumps<1)  return(" Random jumps must be at least one ")   
  }    
   if (!is.null(seed)) {
     if (anyNA(seed) ) return(" NA seed values not allowed") 
     if (any(!is.finite(seed)) | any(is.nan(seed))) return(" Seed must be a finite integer ") 
     if (seed<1)  return(" Seed must be at least one ")   
   } 
  
    rowLevs=rowLevels[(rowLevels>1 | colLevels>1)]
    colLevs=colLevels[(rowLevels>1 | colLevels>1)]
    minplots=prod( replace(rowLevs, rowLevs==1, 2)*replace(colLevs, colLevs==1, 2))
    if ( isTRUE( sum(treatments*replicates) < minplots ) )
      return(paste("The total number of plots is",  sum(treatments*replicates) , 
                   "whereas the required number for the required block design is", minplots))     
    
   if ( isTRUE( sum(treatments) < 2 ) )
     return(paste("The number of treatments must be at least two "))  
   if ( isTRUE( sum(treatments*replicates) < (prod(rowLevels)*2) ) )
     return(paste("There must be at least two plots in each block : the total number of plots is",  sum(treatments*replicates) , 
                  "while the total required number of blocks is", prod(rowLevels),", which is not feasible. "))  
   if ( isTRUE( sum(treatments*replicates) < (prod(rowLevels) + sum(treatments)-1) ) )
     return(paste("The total number of plots is",  sum(treatments*replicates) , 
                  "whereas the total required number of model parameters is", prod(rowLevels) + sum(treatments),", which is not feasible. "))  
   if ( isTRUE( length(rowLevels)!=length(colLevels)) )
     return(paste("The numbers of nested and crossed block level factors must be equal"))  
    
    if ( isTRUE( max(replicates)==2 & min(which((rowLevs==2)==(colLevs==2)))<length(rowLevs)) ) 
      return(paste("For a two replicate treatment design, it is not useful to have nesting within a crossed 2 x 2 row-and-column design "))    

   return(TRUE)
 }
 
 # ******************************************************************************************************************************************************** 
 # Main body of blocks design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
 # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
 # ********************************************************************************************************************************************************     
 if ( is.null(colLevels)) colLevels=rep(1,length(rowLevels))
 testout=testInputs(treatments,replicates,rowLevels,colLevels,searches,seed,jumps) 
 if (!isTRUE(testout)) stop(testout)
 indicator=(rowLevels>1 | colLevels>1)
 rowLevels=rowLevels[indicator]
 colLevels=colLevels[indicator]
 set.seed(seed)
 strata=length(rowLevels)
 # factorial matrices
 z=as.numeric(rbind(rowLevels,colLevels))
 factlevs <- function(r){ gl(z[r],prod(z)/prod(z[1:r])) }
 factMat=do.call(cbind,lapply(1:(2*strata),factlevs))
 hcf=HCF(replicates)
 TF=rep( rep(1:sum(treatments),rep(replicates%/%hcf,treatments)), hcf)
 rand=sample(1:length(TF))
 TF=TF[rand][order(rep(1:hcf,each=length(TF)/hcf)[rand])]
 #block sizes
 blocksizes=sum(treatments*replicates)
 for ( i in 1 :strata) {
  newsizes=NULL
  for (j in 1: length(blocksizes))
    newsizes=c(newsizes,BlockSizes(c(rowLevels[i],colLevels[i]),blocksizes[j] ))
 blocksizes=newsizes 
 }
 # column headers
 colnames=NULL
 for (i in 1 : strata)
   colnames=c(colnames,paste0("Row_",i),paste0("Column_",i))

 Design=data.frame(factMat[rep(1:length(blocksizes),blocksizes),],TF)
 Design[]=lapply(Design, as.factor) 
 colnames(Design)=c(colnames,"Treatments")

 Design=RC_Opt(Design,searches,jumps,strata)
 #add back single rep treatments
 if ( min(replicates)==1 && max(replicates)>1 ) {
   addTF=((sum(treatments[replicates>1])+1) :sum(treatments))
   if (length(addTF)>1) addTF=sample(addTF)
   TF=as.factor(  c(TF, addTF )  ) 
   reptrts=NULL
   for (i in 1 : length(replicates))
     if (replicates[i]>1)
       reptrts=c(reptrts,rep(FALSE,treatments[i])) else 
         reptrts=c(reptrts,rep(TRUE,treatments[i]))
   levels(TF)= (1:sum(treatments*replicates))[order(reptrts)] 
   TF=as.numeric(levels(TF))[TF]
   newblocksizes=Sizes(sum(treatments*replicates),rowLevels)
   Design=data.frame(factMat[rep(1:length(newblocksizes),newblocksizes),])
   TF=TF[order(c(rep(1:length(blocksizes),blocksizes),rep(1:length(blocksizes),(newblocksizes-blocksizes))))]
   blocksizes=newblocksizes
 }
  Design=cbind(Design[,c(1:(2*strata))],rep(1:nrow(Design)),Design[,ncol(Design)])
  Design[]=lapply(Design, as.factor)
  # randomization
  Design=data.frame(do.call(cbind,lapply(1:(2*strata+1), function(r){ sample(nlevels(Design[,r]))[Design[,r]] })) ,Design[,2*strata+2])
  Design=Design[ do.call(order, Design), ]
  Design[]=lapply(Design, as.factor)
  TF=Design[,ncol(Design)]
  
  blocksindex=rep(1,nrow(Design))
  for (i in 1:strata) 
    blocksindex= (blocksindex-1)*nlevels(Design[,2*i-1])*nlevels(Design[,2*i]) + (as.numeric(Design[,2*i-1])-1)*nlevels(Design[,2*i]) + as.numeric(Design[,2*i]) 
   blocksizes=tabulate(blocksindex)
  Design=data.frame(factMat[rep(1:length(blocksizes),blocksizes),] ,TF)
  colnames(Design)=c(colnames,"Treatments")

  if (strata>1) {
    NestStrata=matrix(1,nrow=nrow(Design),ncol=(strata+1))
    for (i in 1:strata) 
      NestStrata[,i+1] =(Design[,2*i-1]-1)*max(Design[,2*i]) + Design[,2*i]
    for (i in 1:strata) 
      NestStrata[,i+1] = (NestStrata[,i]-1)*max(NestStrata[,i+1]) + NestStrata[,i+1]
    for (i in 2:strata) {
      Design[,2*i-1] = (NestStrata[,i]-1)*max(Design[,2*i-1]) + Design[,2*i-1]
      Design[,2*i] =   (NestStrata[,i]-1)*max(Design[,2*i]) + Design[,2*i]  
    }
  }
  colnames(Design)=c(colnames,"Treatments")
  Design[]=lapply(Design, as.factor)
  # aov
  # DF=nestedLevels-1
  # if (strata>1) 
  #   DF[2:strata]=DF[2:strata]*cumblocklevs[1:(strata-1)]
  # if (rowcol & strata>1) DF = c(DF, cumblocklevs[length(cumblocklevs)-1]  *  (min(blocksizes)-1) )
  # if (rowcol & strata==1) DF = c(DF, min(blocksizes)-1 )
  # DF=c(DF,sum(treatments)-1)
  # DF=c(DF,nrow(Design)-sum(DF)-1)
  # AOV=data.frame(DF)
  #  colnames(AOV)="DF"
  #  rownames(AOV)=c(stratumnames,"Treatments","Residual") 
  #print(AOV)
  # efficiencies
  Efficiencies=A_Efficiencies(Design,treatments,replicates,rowLevels,colLevels)
  print(Efficiencies)
  # treatment replications
  Treatments=data.frame(table(Design[,"Treatments"]))
  Treatments[]=lapply(Treatments, as.factor) 
  colnames(Treatments)=c("Treatments","Replicates")
  
    rows=length(blocksizes)
    indicate=c(diff(    as.numeric(Design[,(ncol(Design)-1)])   ,1)==0,FALSE)
    trts=vector( length=rows*cols)
    index=1
    for (i in 1:length(trts)) {
      if (indicate[index]==TRUE) {
        trts[i]=paste(Design[,ncol(Design)][index],Design[,ncol(Design)][index+1])
        index=index+2
      } else {
        trts[i]=Design[,ncol(Design)][index]
        index=index+1
      }
    }
    rc=matrix(trts,nrow=rows,ncol=cols,byrow=TRUE)
    Plan=data.frame(factMat ,rep("",rows),rc)
    colnames(Plan)=c(stratumnames,rep(1:cols)) 
  

    #index=rep(1:length(blocksizes),blocksizes)
    #plots=vector(length=length(blocksizes))
    #l=as.numeric(Design[,"Treatments"])
    #l=formatC(l,width=nchar(max(l)))
    #for (i in 1: length(blocksizes)) 
    #  plots[i]=paste(l[index==i],collapse=" ")
    Plan=data.frame(factMat,rep("",length(blocksizes)) ,plots)
    Plan[is.na(Plan)] = ""
    Plan[]=lapply(Plan,factor) 
    colnames(Plan)=c(stratumnames,"Plots",1:max(blocksizes))
  
  Plan[]=lapply(Plan,as.factor) 
  row.names(Plan)=c(1:nrow(Plan))
  if (strata>1)
  Plan=split(Plan,Plan[1:(strata-1)])
 list(Treatments=Treatments,Efficiencies=Efficiencies,Design=Design,Plan=Plan,AOV=AOV,Seed=seed,Searches=searches,Jumps=jumps) 
} 

