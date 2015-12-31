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
#' The \code{rowblocks} parameter contains the nested blocks for each individual stratum taken in order from the highest to the lowest.
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
#' @param rowblocks factor levels that define the number of nested blocks in each succesive blocks stratum taken in order from the highest to the lowest. 
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
#' blocks(treatments=50,replicates=4,rowblocks=c(4,5))
#' 
#' # as above but with 20 additional single replicate treatments 
#' # giving exactly one single replicate treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),rowblocks=c(4,5))
#' 
#' # 64 treatments x 2 replicates with 2 main blocks and five succesively nested 2-level factors
#' blocks(treatments=64,replicates=2,rowblocks=c(2,2,2,2,2,2))
#' 
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,rowblocks=4)
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
blocks = function( treatments,replicates, rowblocks=HCF(replicates),colblocks=NULL,searches=(1+2000%/%(sum(treatments)+prod(rowblocks))),seed=sample(10000,1),jumps=1) { 
 
 
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
  optEffics=function(TF,BF,ntrts,k) { 
    if (ntrts<=k) 
      e=eigen( (diag(ntrts)-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(ntrts-1)] else    
        e=c(rep(1,(ntrts-k)), 
            eigen((diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])  
    round(c(mean(e)*prod(e/mean(e))^(1/length(e)),1/mean(1/e)),6)
  }
  # ******************************************************************************************************************************************************** 
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
  # ********************************************************************************************************************************************************
  DMax=function(MTT,MBB,MTB,TF,BF,Restrict) {   
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
  Optimise=function(TF,BF,MTT,MBB,MTB,Restrict)  {
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
          if ( !is.na(bound) && isTRUE(all.equal(bound,optEffics(globTF,BF,nlevels(TF),nlevels(BF))[2]))) break
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
        if (counter>1000) return(globTF) # no non-singular swaps
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
  Swaps=function(TF,MF,BF,pivot,rank,n) {
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
    n=length(TF)
    while (rank<fullrank & times<1000) {
      times=times+1
      s=Swaps(TF,MF,BF,pivot,rank,n)
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
    if (times>999)  
      stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design") 
    return(TF)
  }  
  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  rowsOpt=function(TF,MF,BF) { 
    rand=sample(1:length(TF))
    TF=TF[rand][order(MF[rand])]
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
    TF=Optimise(TF,BF,MTT,MBB,MTB,MF)
    TF
  }  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  colsOpt=function(TF,Main,Cols,Rows) { 
    blevels=nlevels(Cols)
    Cols=as.factor((as.numeric(Main)-1)*nlevels(Cols)+as.numeric(Cols))
    BM=Contrasts(Main,Cols)[, rep(c(rep(TRUE,(blevels-1)),FALSE),nlevels(Main)),drop=FALSE]
    TM=Contrasts(Main,TF)[,-nlevels(TF),drop=FALSE] 
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    MBB=matrix(0,nrow=nlevels(Cols),ncol=nlevels(Cols))
    MTT=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))  
    MTB=matrix(0,nrow=nlevels(TF),ncol=nlevels(Cols))
    MBB[1:(nlevels(Cols)-nlevels(Main)), 1:(nlevels(Cols)-nlevels(Main))]=V[1:(nlevels(Cols)-nlevels(Main)),1:(nlevels(Cols)-nlevels(Main)),drop=FALSE]
    MTT[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[(nlevels(Cols)-nlevels(Main)+1):ncol(V),(nlevels(Cols)-nlevels(Main)+1):ncol(V), drop=FALSE]
    MTB[1:(nlevels(TF)-1),  1:(nlevels(Cols)-nlevels(Main)) ]=V[(nlevels(Cols)-nlevels(Main)+1):ncol(V),1:(nlevels(Cols)-nlevels(Main)),drop=FALSE]
    perm=order(order( (1:nlevels(Cols))%%blevels ==0  ))  
    MTB=MTB[,perm]
    MBB=MBB[perm,perm] 
    TF=Optimise(TF,Cols,MTT,MBB,MTB,Rows)
    TF
  }  
  
  # ******************************************************************************************************************************************************** 
  # Generates an initial orthogonal design then builds algebraic lattice blocks or calls the general block design algorithm as appropriate
  # ********************************************************************************************************************************************************     
  initializeTF=function(treatlevs,replevs) {
    
    
    # special balanced designs 
    
    v=sqrt(ntrts)  
      regrep=isTRUE(all.equal(max(replevs),min(replevs)))
      
      for (stratum in 1 : length(rowblocks)) { 
        nblocks=prod(rowblocks[1:stratum])
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
        } else TF=rowsOpt(TF,Design,stratum,rowblocks,hcf)
      }
     
    levels(TF)=sample(1:ntrts) 
    TF
  }
  
  
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for the blocks in each stratum of a design 
  # ********************************************************************************************************************************************************     
  A_Efficiencies=function(Design,treatments,replicates,ncols) {
    strata=ncol(Design)-2
    Columns= Design[,strata+1]
    if (strata>1) Columns = as.factor((as.numeric(Design[,strata-1])-1)*nlevels(Columns) +  as.numeric(Columns))
    hcf=HCF(replicates)
    nunits=nrow(Design)
    ntrts=nlevels(Design[,strata+2])
    nblocks=as.numeric(sapply(Design,nlevels)[1:strata])
    if (ncols>1) nblocks=c(nblocks,nlevels(Columns))
    if (ncols==1) effics=matrix(1,nrow=strata,ncol=2) else effics=matrix(1,nrow=(strata+1),ncol=2)
    if (ncols==1) bounds=rep(NA,strata) else bounds=rep(NA,(strata+1))
    for (i in 1:strata) { 
      if (isTRUE(all.equal(max(replicates),min(replicates))) ) {
        if ( isTRUE(all.equal(nunits%%nblocks[i],0))) 
          bounds[i]=upper_bounds(nunits,sum(treatments),nblocks[i]) else if (hcf%%nblocks[i]==0) 
            bounds[i]=1
      }
      if ( sum(treatments)>1 && nblocks[i]>1)
        effics[i,]=optEffics(Design$Treatments,Design[,i],ntrts,nblocks[i])  
    }
    if (ncols>1){
      bounds[strata+1]=upper_bounds(nunits,sum(treatments),nblocks[strata+1]) 
      effics[strata+1,]=optEffics(Design$Treatments,Columns,ntrts,nblocks[strata+1]) 
    }
    names=names(Design)[1:strata]
    if (ncols>1) names=c(names,"Columns")
    efficiencies=data.frame(cbind(names,nblocks, effics, bounds)) 
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies[,'Blocks'] = as.factor(efficiencies[,'Blocks'])
    efficiencies
  }
  # ******************************************************************************************************************************************************** 
  # Carries out some input validation
  # ********************************************************************************************************************************************************     
  testInputs=function(treatments,replicates,rowblocks,seed) {  
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
    
    if (!is.null(rowblocks)) {
      if (anyNA(rowblocks) ) return(" NA rowblocks values not allowed") 
      if (any(!is.finite(rowblocks)) | any(is.nan(rowblocks)) ) return(" rowblocks can contain only finite integers ")
      if (min(rowblocks)<1) return (" rowblocks must be at least one ")
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
    if ( isTRUE( sum(treatments) < 2 ) )
      return(paste("The number of treatments must be at least two "))  
    
    if ( isTRUE( sum(treatments*replicates) < (prod(rowblocks) + sum(treatments)-1) ) )
      return(paste("The total number of plots is",  sum(treatments*replicates) , 
                   "whereas the total required number of model parameters is", prod(rowblocks) + sum(treatments),", which is not feasible. "))  
    return(TRUE)
  }
  
  # ******************************************************************************************************************************************************** 
  # Main body of blocks design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************     

  testout=testInputs(treatments,replicates,rowblocks,seed) 
  if (!isTRUE(testout)) stop(testout)
  set.seed(seed)
  if (isTRUE(all.equal(length(rowblocks),0))) rowblocks=1
  if (isTRUE(all.equal(length(colblocks),0))) colblocks=rep(1,length(rowblocks))
  if (!isTRUE(all.equal(length(colblocks),length(rowblocks)))) stop("The number of row blocks strata and the number of column blocks strata must be equal")
  fulltreatments=treatments
  fullreplicates=replicates
  treatments=treatments[replicates>1]
  replicates=replicates[replicates>1]
  hcf=HCF(replicates)
  MB=rep(1:hcf,each=(sum(treatments*replicates)/hcf) )
  TF=rep( rep(1:sum(treatments),rep(replicates/hcf,treatments)), hcf)
  rand=sample(sum(treatments*replicates))
  TF=as.factor( TF[rand][order(MB[rand])] )
  #colblocks=colblocks[rowblocks>1]
  #rowblocks=rowblocks[rowblocks>1]
  strata=length(rowblocks)
  cumrowlevs=cumprod(rowblocks)
  cumcollevs=cumprod(colblocks)

  stratumnames=sapply(1:strata, function(i) { paste("Rows",i,"*","Columns", i) })
  
  

  rowfacMat= matrix(nrow=prod(rowblocks[1:strata]),ncol=strata)
  rowfactlevs = function(r){ gl(cumrowlevs[r],prod(rowblocks[1:strata])/cumrowlevs[r]) }
  rowfacMat=do.call(cbind,lapply(1:strata,rowfactlevs))
  
  colfacMat= matrix(nrow=prod(colblocks[1:strata]),ncol=strata)
  colfactlevs = function(r){ gl(cumcollevs[r],prod(colblocks[1:strata])/cumcollevs[r]) }
  colfacMat=do.call(cbind,lapply(1:strata,colfactlevs))
  

  blocks=rep(1,blocksizes)
  strata=length(rowblocks)
  cumblocks=c(1,cumprod(rowblocks*colblocks))
  design=matrix(nrow=blocksizes,ncol=2*strata)
  for (i in 1 :strata) {
    v=vector(mode = "list", length =length(blocksizes))
    for (j in 1:length(blocksizes)) { 
      v[[j]]=rep(floor(blocksizes[j]/rowblocks[i]),rowblocks[i])
      resid=blocksizes[j]-sum(v[[j]])
      if (resid>0)
        v[[j]][1:resid]=v[[j]][1:resid]+1
    }
    rowsizes=unlist(v)
    rows=rep(1:length(rowsizes),rowsizes)
    w=vector(mode = "list", length =length(rowsizes))
    shift=0
    for (j in 1:length(rowsizes)) {
      colsbase=floor(rowsizes[j]/colblocks[i])
      resid=rowsizes[j]-colblocks[i]*colsbase
      w[[j]]=rep(colsbase,colblocks[i])
      if (resid>0) {
        w[[j]][(shift:(shift+resid-1))%%colblocks[i]+1]=colsbase+1
        shift=shift+resid
      } 
    }
    blocksizes=unlist(w)
    rowfactlevs = rep(rep(1:rowblocks[i],each=colblocks[i]),cumblocks[i])
    colfactlevs = rep(rep(1:colblocks[i],rowblocks[i]),cumblocks[i])
    blocklevels=  rep(1:cumblocks[i],each=colblocks[i]*rowblocks[i])
    rowfactlevs=(blocklevels-1)*rowblocks[i] + rowfactlevs
    colfactlevs=(blocklevels-1)*colblocks[i] + colfactlevs
    # print(cbind(blocklevels,rowfactlevs,colfactlevs,blocksizes))
  }
  
  
  Design=as.data.frame(design)
  Design[]=lapply(Design, as.factor) 
  

  for ( i in 1: strata) {
 
    
    TF=rowsOpt(TF,Main,Rows)
    TF=colsOpt(TF,Main,Cols,Rows)
  }
    
    
 



  #  Design[,ncol(Design)-1]=as.factor(unlist(v))
   # colnames(Design)=c(stratumnames,"Columns","Treatments") 
    #optimize columns design
   # Design=colsOpt(Design)
  #  TF=colsOpt(TF,Main,Cols,Rows)
 # }
  
  
  
  #plan=sapply(v, tabulate)
  #print(plan)
  # efficiencies
  Efficiencies=A_Efficiencies(Design,treatments,replicates,colblocks[strata])
  
  # aov
  DF=rowblocks-1
  if (strata>1) 
    DF[2:strata]=DF[2:strata]*cumrowlevs[1:(strata-1)]
  DF=c(DF,sum(treatments)-1)
  DF=c(DF,nrow(Design)-sum(DF)-1)
  AOV=data.frame(DF)
  colnames(AOV)="DF"
  rownames(AOV)=c(stratumnames,"Treatments","Residual")
  
  # treatment replications
  Treatments=data.frame(table(Design[,"Treatments"]))
  Treatments[]=lapply(Treatments, as.factor) 
  colnames(Treatments)=c("Treatments","Replicates")
 
  # Factorial block levels
  Design[,c(1:strata)] = do.call(cbind,lapply(1:strata, function(r){ (as.numeric(Design[,r])-1)%%rowblocks[r]+1 }))
  Design[]=lapply(Design, as.factor)

 # Plan layout with rows or blocks arranged vertically and plots or columns arranged horizontally 
  #print(v)
  
  if (colblocks[strata]>1) {
    rc = matrix("",nrow=cumrowlevs[strata],ncol=colblocks[strata],byrow=TRUE)
    plot = 1
    for ( i in 1:nrow(rc)) {
      for (icol in v[[i]]) {
        rc[i,icol]=paste(rc[i,icol],formatC(as.numeric( Design[plot,ncol(Design)],width=nchar(sum(treatments)))))
        plot=plot+1
      }
  }
}
  if (colblocks[strata]==1) {
    if (min(rowsizes)<max(rowsizes)) {
      rc = matrix(formatC("-",width=nchar(sum(treatments))),nrow=cumrowlevs[strata],ncol=max(rowsizes),byrow=TRUE)
      plot = 1
      for ( i in 1:nrow(rc)) {
        for (j in 1:rowsizes[i]) {
          rc[i,j]=formatC(as.numeric(Design[plot,ncol(Design)]),width=nchar(sum(treatments)))
          plot=plot+1
        }
      }
    } else 
      rc=matrix( formatC(as.numeric(Design[,ncol(Design)] ,width=nchar(sum(treatments)))),nrow=cumrowlevs[strata],ncol=min(rowsizes),byrow=TRUE)
  }
  rowfacMat[,c(1:strata)] = do.call(cbind,lapply(1:strata, function(r){ (as.numeric(rowfacMat[,r])-1)%%rowblocks[r]+1 }))
  Plan=data.frame( cbind(rowfacMat, rep("",nrow(rc)), rc))
  if (colblocks[strata]>1) colnames(Plan)=c(stratumnames,"Columns",1:ncol(rc))
  if (colblocks[strata]==1) colnames(Plan)=c(stratumnames,"Plots",1:ncol(rc))
  Plan[]=lapply(Plan,as.factor) 
  if (strata>1 & colblocks[strata]>1)
   Plan=split(Plan,Plan[1:(strata-1)]) 
  else if (strata>1 & colblocks[strata]==1)
    Plan=split(Plan,Plan[1:(strata)])
  
  list(Treatments=Treatments,Efficiencies=Efficiencies,Design=Design,Plan=Plan,AOV=AOV,Seed=seed,Searches=searches,Jumps=jumps) 
} 