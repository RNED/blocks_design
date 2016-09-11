#' @title Block designs
#'
#' @description
#'
#' Constructs randomized nested and crossed block designs for unstructured treatment sets where treatments can have any arbitrary levels of replication
#' and blocks can have any arbitrary feasible depth of nesting.
#'
#' @details
#'
#' The \code{blocksdesign} package constructs arbitrary block designs with arbitrary depth of nesting and arbitrary crossed row-and-column blocks
#' in each nested stratum. Where one classification of a row-and-column design has a single block level, the design automatically defaults to a simple nested blocks design
#' (if both classifications of a stratum have a single level that stratum is dropped from the design).
#'
#'  The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order. \cr
#'  \item  A table showing the replication number of each treatment in the design. \cr
#'  \item  A table showing the block levels and the achieved D- and A-efficiency factors for each blocks stratum together with A-efficiency upper bounds, where available. \cr
#'  \item  Plans showing the allocation of treatments to blocks or to rows and to columns in the bottom stratum of the design. \cr
#' }
#'
#' @param treatments  a vector giving a partition of the total required number of treatments into sets of equally replicated treatments.
#'
#' @param replicates  a vector giving the replication number of each equally replicated treatment set in the \code{treatments} vector.
#'
#' @param rows  a vector of factor levels for the row blocks in each succesive stratum of the blocks design taken in order from the highest to the lowest.
#' The default is a single set of main blocks equal to the hcf of the replication numbers.
#'
#' @param columns  a vector of factor levels for the column blocks in each succesive stratum of the blocks design taken in order from the highest to the lowest.
#' The \code{rows} and the \code{columns} vectors, if both present, must be of equal length. The default is the null vector.
#'
#' @param seed  an integer initializing the random number generator. The default is a random seed.
#'
#' @param searches  the maximum number of local optima searched for a design optimization. The default is 1 plus the floor of 10000 divided by the number of plots.
#'
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#'
#' @return
#' \item{Design}{Data frame giving the optimized block and treatment factors in plot order}
#' \item{Plan}{Data frame giving a plan view of the treatments design in the bottom stratum of the design classified by rows and columns}
#' \item{Efficiencies}{The achieved A- and D-efficiencies for each stratum of the design together with an A-efficiency upper-bound, where available}
#' \item{seed}{Numerical seed for random number generator}
#' \item{searches}{Maximum number of searches in each stratum}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima}
#'
#' @references
#'
#' @examples
#'
#' @export
#' @importFrom stats anova lm
#'
#'
factblocks = function(TF,replicates=1,rows=NULL,columns=NULL,order=NULL,searches=(1+10000%/%nrow(TM)),seed=sample(10000,1),jumps=1) {
  options(contrasts=c('contr.sum','contr.poly'))
  if (!is.data.frame(TF)) stop("Treatment factors not in a valid data frame")
  fInd=sapply(TF, is.factor)
  fTF=TF[, fInd,drop=FALSE]
  pTF=TF[,!fInd,drop=FALSE]

  if (ncol(fTF)>0) 
    for (i in 1: min(order,ncol(fTF))) {
      comb=combn(ncol(fTF),i)
      pTF=c(pTF,prod(pTF[,comb]))
    }
  
  #lapply(poly(pTF,degree=4))
  P=lapply(1:ncol(pTF), function(r){ poly(pTF[,r],degree=min(order,length(unique(pTF[,r]))-1),raw=TRUE)})

  sapply(fTF, nlevels)
  length(unique(pTF[[]]))
  
  model.matrix(~ a * poly(TFP,degree=4) , TF)
  TM=model.matrix(  as.formula(paste(" ~ ", paste(colnames(TF), collapse= "*")))  ,TF)
print(TM)

  # ********************************************************************************************************************************************************
  # Finds the highest common factor (hcf) of a set of numbers omitting any zero values (Euclidean algorithm)
  # ********************************************************************************************************************************************************
  HCF=function(replicates)  {
    replicates=sort(replicates[replicates>0])
    for (i in  seq_len(length(replicates)))
      while (!isTRUE(all.equal(replicates[i]%%replicates[1],0))) replicates[c(1,i)] = c(replicates[i]%%replicates[1], replicates[1])
      replicates[1]
  }
  # ********************************************************************************************************************************************************
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
  # ********************************************************************************************************************************************************
  UpDate=function(MTT,MBB,MTB,t,b) {
    t=c(t,0)
    b=c(b,0)
    Mttt=crossprod(MTT,t)
    Mbbb=crossprod(MBB,b)
    Mtbt=crossprod(MTB,t)
    Mtbb=crossprod(t(MTB),b)
    tMt=as.numeric(crossprod(t,Mttt))
    bMb=as.numeric(crossprod(b,Mbbb))
    tMb=as.numeric(crossprod(b,Mtbt))
    f1=(Mttt+Mtbb)/sqrt(2)
    f2=(Mbbb+Mtbt)/sqrt(2)
    g1=(Mtbb-Mttt)/sqrt(2)
    g2=(Mbbb-Mtbt)/sqrt(2)
    a=(tMt+bMb+2*tMb)/2
    b=(tMt+bMb-2*tMb)/2
    c=(bMb-tMt)/2
    d=(1+a)*(1-b)+c*c
    MTT=MTT- (tcrossprod(f1)*(1-b) - tcrossprod(g1)*(1+a) + (tcrossprod(g1,f1)+tcrossprod(f1,g1))*c)/d
    MBB=MBB- (tcrossprod(f2)*(1-b) - tcrossprod(g2)*(1+a) + (tcrossprod(g2,f2)+tcrossprod(f2,g2))*c)/d
    MTB=MTB- (tcrossprod(f1,f2)*(1-b) - tcrossprod(g1,g2)*(1+a) + (tcrossprod(g1,f2)+tcrossprod(f1,g2))*c)/d
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }


  # ********************************************************************************************************************************************************
  # Random swaps
  # ********************************************************************************************************************************************************
  Swaps=function(TM,MF,BF,pivot,rank,CrossedRestrict) {
    n=nrow(TM)
    candidates=NULL
    while (isTRUE(all.equal(length(candidates),0))) {
      if (rank<(n-1)) s1=sample(pivot[(1+rank):n],1) else s1=pivot[n]
      candidates = seq_len(n)[MF==MF[s1] & CrossedRestrict==CrossedRestrict[s1] & BF!=BF[s1] & TM!=TM[s1,]]
    }
    if ( length(candidates)>1 )
      s2=sample(candidates,1) else s2=candidates[1]
      s=c(s1,s2)
  }
  # ********************************************************************************************************************************************************
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************
  NonSingular=function(TM,MF,BF,CrossedRestrict) {
    BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
    BM[cbind(1:length(BF),BF)]=1
    fullrank=ncol(TM)+ncol(BM)
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    times=0
    while (rank<fullrank & times<1000) {
      times=times+1
      s=Swaps(TM,MF,BF,pivot,rank,CrossedRestrict)
      rindex=seq_len(nrow(TM))
      rindex[c(s[1],s[2])]=rindex[c(s[2],s[1])]
      newQ=qr(t(cbind(BM,TM[rindex,])))
      if (isTRUE(all.equal(newQ$rank,rank)) || newQ$rank>rank) {
        TM=TM[rindex,]
        rank=newQ$rank
        pivot=newQ$pivot
      }
    }
    if (times>999) TM=NULL
    return(TM)
  }
  # ********************************************************************************************************************************************************
  # Contrasts for factor NF centered within the levels of factor MF to ensure that NF information is estimated within the levels of factor MF only
  # ********************************************************************************************************************************************************
  Contrasts=function(MF,NF) {
    NM=matrix(0,nrow=length(NF),ncol=nlevels(NF))
    NM[cbind(seq_len(length(NF)),NF)]=1 # factor indicator matrix
    do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(NM[MF==i,] , center = TRUE, scale = FALSE)}))
  }
  # ********************************************************************************************************************************************************
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
  # ********************************************************************************************************************************************************
  DMax=function(MTT,MBB,MTB,TF,Restrict,BF,Blocks,TM,BM){
    relD=1
    mainSizes=tabulate(Restrict)
    nSamp=pmin( rep(8,nlevels(Restrict)), mainSizes)
    repeat {
      improved=FALSE
      for (k in 1: nlevels(Restrict)) {
        s=sort(sample( seq_len(length(TF)) [Restrict==k], nSamp[k]))
        TB=MTB[TF[s],BF[s],drop=FALSE]-tcrossprod(MTB[cbind(TF[s],BF[s])],rep(1,nSamp[k]))
        dMat=(TB+t(TB)+1)**2-
          (2*MTT[TF[s],TF[s],drop=FALSE]- tcrossprod(MTT[cbind(TF[s],TF[s])]+rep(1,nSamp[k])) + tcrossprod(MTT[cbind(TF[s],TF[s])]) + 1)*
          (2*MBB[BF[s],BF[s],drop=FALSE]- tcrossprod(MBB[cbind(BF[s],BF[s])]+rep(1,nSamp[k])) + tcrossprod(MBB[cbind(BF[s],BF[s])]) + 1)
        is.na(dMat[ dMat<0.00001|is.na(dMat) ] ) = TRUE
        sampn=which.max(dMat)
        i=1+(sampn-1)%%nSamp[k]
        j=1+(sampn-1)%/%nSamp[k]
        maxd=dMat[i,j]
        if  (  maxd>1  && dMat[i,j]>1  &&  !isTRUE(all.equal(dMat[i,j],1))  ) {
          improved=TRUE
          relD=relD*maxd
          t=TM[s[i],]-TM[s[j],]
          b=BM[s[j],]-BM[s[i],]
          up=UpDate(MTT,MBB,MTB,t,b)
          MTT=up$MTT
          MBB=up$MBB
          MTB=up$MTB
          TM[c(s[i],s[j]),]=TM[c(s[j],s[i]),]
          TF[c(s[i],s[j])]=TF[c(s[j],s[i])]
        }
      }
      if (improved) next
      if (sum(nSamp) < min(length(TF),512)) nSamp=pmin(mainSizes,2*nSamp) else break
    }
    list(MTT=MTT,MBB=MBB,MTB=MTB,TF=TF,TM=TM,relD=relD)
  }
  # ********************************************************************************************************************************************************
  #  Searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # ********************************************************************************************************************************************************
  Optimise=function(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,TM,BM) {
    if (is.null(Rows)) {
      Restrict=rep(1,length(TF))
      BF=Main
    } else if (is.null(Columns)) {
      Restrict=Main
      BF=Rows
    } else {
      Restrict=Rows
      BF=Columns
    }
    globrelD=0
    relD=1
    globTF=TF
    breps=tabulate(BF)
    for (r in 1:searches) {
      dmax =  DMax(MTT,MBB,MTB,TF,Restrict,BF,Blocks,TM,BM)
      if ( !isTRUE(all.equal(dmax$relD,1)) && dmax$relD>1) {
        relD=relD*dmax$relD
        TF=dmax$TF
        TM=dmax$TM
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB
        if (!isTRUE(all.equal(relD,globrelD)) && relD>globrelD) {
          globTF=TF
          globrelD=relD
        }
      }
      if (r==searches) break
      for (iswap in 1:jumps) {
        counter=0
        repeat {
          counter=counter+1
          s1=sample(seq_len(length(TF)),1)
          z=seq_len(length(TF))[Restrict==Restrict[s1] & BF!=BF[s1] & TF!=TF[s1]]
          if (length(z)==0) next
          if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z[1])
          Dswap = (1+MTB[TF[s[1]],BF[s[2]]]+MTB[TF[s[2]],BF[s[1]]]-MTB[TF[s[1]],BF[s[1]]]-MTB[TF[s[2]],BF[s[2]]])**2-
            (2*MTT[TF[s[1]],TF[s[2]]]-MTT[TF[s[1]],TF[s[1]]]-MTT[TF[s[2]],TF[s[2]]])*(2*MBB[BF[s[1]],BF[s[2]]]-MBB[BF[s[1]],BF[s[1]]]-MBB[BF[s[2]],BF[s[2]]])
          if (Dswap>.1 | counter>1000) break
        }
        if (counter>1000) return(globTF) # no non-singular swaps
        relD=relD*Dswap
        t=TM[s[1],]-TM[s[2],]
        b=BM[s[2],]-BM[s[1],]
        up=UpDate(MTT,MBB,MTB,t,b)
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]
        TF[c(s[1],s[2])]=TF[c(s[2],s[1])]
      }
    }
    globTF
  }
  # ******************************************************************************************************************************************************** 
  # Optimize the nested columns blocks assuming a possible set of Main block constraints Initial randomized starting design. 
  # If the row-column intersections contain 2 or more plots a weighted (columns + w*rows.columns) model is fitted for w>=0 and w<1 
  # ********************************************************************************************************************************************************    
  colsOpt=function(TF,Main,Rows,Columns,weighted) { 
    TF=NonSingular(TF,Main,Columns,Rows)
    if (is.null(TF)) return(TF)
    main=nlevels(Main)
    ncol=nlevels(Columns)/main
    nrow=nlevels(Rows)/main
    nblocks=nrow*ncol
    blocks=main*nblocks
    CM=Contrasts(Main,Columns)[,rep(c(rep(TRUE,(ncol-1)),FALSE),main),drop=FALSE]
    TM=Contrasts(Main,TF)[,-ntrts,drop=FALSE] 
    V=chol2inv(chol(crossprod(cbind(CM,TM))))
    indicv=seq(ncol(CM)+1,ncol(TM)+ncol(CM))
    MTT=rbind(cbind(V[indicv,indicv,drop=FALSE],rep(0,(ntrts-1))),rep(0,ntrts))
    MBB=matrix(0,nrow=nlevels(Columns),ncol=nlevels(Columns))
    MTB=matrix(0,nrow=ntrts,ncol=nlevels(Columns))
    MBB[seq_len(ncol(CM)),seq_len(ncol(CM))]=V[seq_len(ncol(CM)),seq_len(ncol(CM)),drop=FALSE]
    MTB[seq_len(ncol(TM)),seq_len(ncol(CM))]=V[indicv,seq_len(ncol(CM)),drop=FALSE]
    perm=c(rbind(matrix(seq_len(ncol(CM)),nrow=ncol-1,ncol=main),seq_len(main)+ncol(CM)))
    MTB=MTB[,perm]
    MBB=MBB[perm,perm] 
    Blocks=as.factor((as.numeric(Main)-1)*nblocks + ((as.numeric(Rows)-1)%%nrow)*ncol + (as.numeric(Columns)-1)%%ncol + 1)
    BM=Contrasts(Main,Blocks)[,rep( c(rep(TRUE,(nblocks-1)),FALSE),main),drop=FALSE]
    DM=cbind(BM,TM)
    if ( length(TF)>=(blocks+nlevels(TF)) && qr(t(DM))$rank==ncol(DM) && isTRUE(weighted))  {
      V=chol2inv(chol(crossprod(DM)))
      indicv=seq_len(ncol(TM))+ncol(BM)
      Mtt=rbind(cbind(V[indicv,indicv,drop=FALSE],rep(0,(ntrts-1))),rep(0,ntrts))
      Mbb=matrix(0,nrow=blocks,ncol=blocks)
      Mtb=matrix(0,nrow=ntrts,ncol=blocks)
      Mbb[seq_len(ncol(BM)),seq_len(ncol(BM))]=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
      Mtb[seq_len(ncol(TM)),seq_len(ncol(BM))]=V[indicv,seq_len(ncol(BM)),drop=FALSE]
      perm=c(rbind(matrix(seq_len(ncol(BM)),nrow=nblocks-1,ncol=main),seq_len(main)+ncol(BM)))
      Mtb=Mtb[,perm]
      Mbb=Mbb[perm,perm] 
    } else {
      Mtb=NULL
      Mtt=NULL
      Mbb=NULL
    }
    TF=Optimise(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,Mtt,Mbb,Mtb,weighted)
    return(TF)
  }  
  # *******************************************************************************************************************************************************
  # Optimize the nested Blocks assuming a possible set of Main block constraints Initial randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************
  rowsOpt=function(TM,Main,Rows) {
    main=nlevels(Main) # main blocks
    blocks=nlevels(Rows)
    TM=NonSingular(TM,Main,Rows,rep(1,nrow(TM)))
    if (is.null(TM)) return(TM)
    BM=Contrasts(Main,Rows)[, rep(c(rep(TRUE,((blocks/main)-1)),FALSE),main),drop=FALSE]
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    indicv=seq_len(ncol(TM))+ncol(BM)
    MTT=rbind(cbind(V[indicv,indicv,drop=FALSE],rep(0,(ntrts-1))),rep(0,ntrts))
      MBB=matrix(0,nrow=blocks,ncol=blocks)
      MTB=matrix(0,nrow=ntrts,ncol=blocks)
      MBB[seq_len(ncol(BM)),seq_len(ncol(BM))]=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
      MTB[seq_len(ncol(TM)),seq_len(ncol(BM))]=V[indicv,seq_len(ncol(BM)),drop=FALSE]
      perm=c(rbind(matrix(seq_len(ncol(BM)),nrow=blocks/main-1,ncol=main),seq_len(main)+ncol(BM)))
      MTB=MTB[,perm]
      MBB=MBB[perm,perm]
      TF=Optimise(TF,Main,Rows,NULL,NULL,MTT,MBB,MTB,TM,BM)
    return(TF)
  }
  # ********************************************************************************************************************************************************
  # Finds row and column sizes in each stratum of a design
  # ********************************************************************************************************************************************************
  Sizes=function(blocksizes,stratum) {
    nblocks=length(blocksizes)
    newblocksizes=NULL
    for (j in seq_len(nblocks)) {
      rowsizes=rep(blocksizes[j]%/%rows[stratum],rows[stratum])
      resid=blocksizes[j]-sum(rowsizes)
      if (resid>0)
        rowsizes[1:resid]=rowsizes[1:resid]+1
      rowcolsizes=vector(mode = "list", length =rows[stratum])
      for ( z in 1:rows[stratum])
        rowcolsizes[[z]]=rep(rowsizes[z]%/%columns[stratum] , columns[stratum])
      shift=0
      for (z in seq_len(rows[stratum])) {
        resid=rowsizes[z]-sum(rowcolsizes[[z]])
        if (resid>0) {
          rowcolsizes[[z]][(shift:(shift+resid-1))%%columns[stratum]+1]=rowcolsizes[[z]][(shift:(shift+resid-1))%%columns[stratum]+1]+1
          shift=shift+resid
        }
      }
      newblocksizes=c(newblocksizes,unlist(rowcolsizes))
    }
    newblocksizes
  }
  # ********************************************************************************************************************************************************
  # Main body of rows design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************
  TM=TM[,-1]
  TM=do.call("rbind", rep(list(TM), replicates))
  if (length(columns)==0) columns=rep(1,length(rows))
  indic=rows*columns
  if (max(indic)==1) {
    rows=1
    columns=1
  } else {
    rows=rows[(indic!=1)]
    columns=columns[(indic!=1)]
  }
  cumcols=cumprod(columns)
  cumrows=cumprod(rows)
  cumblocks=c(1,cumprod(rows*columns))
  strata=length(rows)
  plots=nrow(TM)
  set.seed(seed)
    # Design factors
    rowcol=as.numeric(rbind(rows,columns))
    cpblocks=c(1,cumprod(rows*columns))
    cprowcol=cumprod(rowcol)
    isrowcol=max(columns>1)
    fDesign=do.call(cbind,lapply(1:(2*strata),function(i) {  gl(  rowcol[i],   cprowcol[2*strata]/cprowcol[i], cprowcol[2*strata]  )    })) -1
    BlocksInStrata=do.call(cbind,lapply(1:(strata+1),function(i) {  gl(  cpblocks[i], cpblocks[strata+1]/cpblocks[i], cpblocks[strata+1]  )    })) -1
    if (isrowcol)  fDesign=data.frame(do.call(cbind,lapply(1:(2*strata), function(r){fDesign[,r]+BlocksInStrata[,(r-1)%/%2+1]*rowcol[r] }))+1)
    if (!isrowcol) fDesign=data.frame(do.call(cbind,lapply(1:strata,     function(r){fDesign[,2*r-1]+BlocksInStrata[,r]*rowcol[2*r-1] }))+1)
    fDesign[]=lapply(fDesign, as.factor)
    BlocksInStrata=data.frame(BlocksInStrata+1)
    BlocksInStrata[]=lapply(BlocksInStrata, as.factor)
    #permBlocks will randomise whole blocks in nested strata or will randomize rows and columns in each nested stratum of a crossed design
    tempDesign=data.frame( do.call(cbind,lapply(1:ncol(fDesign), function(r){ sample(nlevels(fDesign[,r]))[fDesign[,r]] })) , seq_len(nrow(fDesign)   ))
    permBlocks=tempDesign[ do.call(order, tempDesign), ][,ncol(tempDesign)]
    nunits=nrow(TM)
    # design
    regular=TRUE
    blocksizes=nunits
    for (i in seq_len(strata)) {
      blocksizes=Sizes(blocksizes,i)
      if (max(blocksizes) > min(blocksizes)) regular=FALSE
    }
     Design  = data.frame(fDesign[rep(seq_len(nrow(fDesign)),  blocksizes ),])
    BlocksInStrata  = data.frame(BlocksInStrata[rep(seq_len(nrow(BlocksInStrata)),  blocksizes ),])
    Design[]= lapply(Design, as.factor)
    BlocksInStrata[]= lapply(BlocksInStrata, as.factor)
    hcf=HCF(replicates)
    print(TM)
    for ( z in seq_len(10)) {
      rand=sample(nunits)
      TM=TM[rand,][order(rep(seq_len(hcf),each=nunits/hcf)[rand]),] 
      for ( i in seq_len(strata)) {
        if (!isrowcol && rows[i]>1 && !all(hcf%%cumprod(rows[1:i])==0)) TM=rowsOpt(TM,BlocksInStrata[,i],Design[,i])
        if (is.null(TM)) break
        if (isrowcol &&  rows[i]>1 && !all(hcf%%cumprod(rows[1:i])==0)) TM=rowsOpt(TM,BlocksInStrata[,i],Design[,2*i-1])
        if (is.null(TM)) break
        if (isrowcol && columns[i]>1) TM=colsOpt(TM,BlocksInStrata[,i],Design[,2*i-1],Design[,2*i])
        if (is.null(TM)) break
      }
      if (!is.null(TM)) break
    }
print(TM)

    if (is.null(TM)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
    # add back single rep treatments for nested stratum blocks only
    
    if ( min(replicates)==1 && max(replicates)>1 && ( max(columns)==1 || regular==TRUE ) ) {
      nunits=nrow(TM)
      ntrts=sum(treatments)
      fullblocksizes=nunits
      for (i in seq_len(strata))
        fullblocksizes=Sizes(fullblocksizes,i)
      sblocksizes=fullblocksizes-blocksizes
      TF=as.numeric(levels(TF))[TF]
      addTF=rep(1:sum(treatments))[fullreps==1]
      BlocksInStrata=rep(1:length(blocksizes),blocksizes)
      repTF=split(TF,BlocksInStrata)
      addBlocks=which(sblocksizes>0)
      sblocks=rep( 1:length(sblocksizes),sblocksizes)
      singTF=split(addTF,sblocks)
      for (i in 1:length(addBlocks))
        repTF[[addBlocks[i]]]=sample(append( repTF[[addBlocks[i]]] ,  singTF[[i]]))
      TF=as.factor(unlist(repTF))
      blocksizes=fullblocksizes
    }
    # Randomize

    D=as.data.frame(cbind(rep(1:length(blocksizes),blocksizes),sample(seq_len(nunits)),TM))
    D[]=lapply(D, as.factor)
    D[,1] = factor(D[,1],levels(D[,1])[permBlocks])
    TF=D[ do.call(order, D),][,ncol(D)]
    blocksizes=blocksizes[permBlocks]
    if (isrowcol)
      stratumnames=unlist(lapply(1:strata, function(i) { c(paste("Rows",i), paste("Columns", i) )})) else
        stratumnames=unlist(lapply(1:strata, function(i) {paste("Blocks",i)}))
    colnames(fDesign)=stratumnames
    Design =fDesign[rep(seq_len(nrow(fDesign)),  blocksizes ),]
    Design =as.data.frame(cbind( Design ,TF) )
    Design[]  = lapply(Design, as.factor)
    colnames(Design)=c(stratumnames,"Treatments")
    rownames(Design)=NULL
    #Plan
    V=split(Design[,ncol(Design)],rep(1:cumblocks[strata+1],blocksizes))
    if (!isrowcol| columns[length(columns)]==1  | rows[length(rows)]==1  ) {
      Plots=rep("",length(V))
      Plan=as.data.frame(cbind(fDesign,Plots , do.call(rbind, lapply(V, function(x){ length(x) =max(blocksizes); x }))))
    } else {
      if(strata>1) {
        fDesign=do.call(cbind,lapply(1:(2*(strata-1)),function(i) {gl(rowcol[i],cprowcol[2*(strata-1)]/cprowcol[i], cprowcol[2*(strata-1)])}))
        fDesign= data.frame(cbind( fDesign[ rep(seq(nrow(fDesign)),each=rows[strata]), ], seq_len(nrow(fDesign))))
      } else fDesign=data.frame(seq_len(rows[1]))
      colnames(fDesign)=stratumnames[1:ncol(fDesign)]
      rownames(fDesign)=NULL
      fDesign[]=lapply(fDesign, as.factor)
      rcblocks=unlist(lapply(V, paste, collapse = ",") )
      plan=matrix("",nrow=cumblocks[strata]*columns[strata],ncol=cumblocks[strata]*rows[strata])
      for (z in 1: cumblocks[strata])
        plan[c(((z-1)*columns[strata]+1):(z*columns[strata])),c(((z-1)*rows[strata]+1):(z*rows[strata]))] =
        rcblocks[c(((z-1)*rows[strata]*columns[strata]+1):(z*rows[strata]*columns[strata]))]
      Columns=rep("",nrow(fDesign))
      Plan=as.data.frame(cbind(fDesign,Columns,t(plan)))
      names(Plan)[names(Plan) == 'Columns'] = paste('Columns', length(columns))
    }
  # treatment replications
  TreatmentsTable=data.frame(table(Design[,"Treatments"]))
  TreatmentsTable[]=lapply(TreatmentsTable, as.factor)
  colnames(TreatmentsTable)=c("Treatments","Replicates")
  list(Treatments=TreatmentsTable,Plan=Plan,Design=Design,Seed=seed,Searches=searches,Jumps=jumps)
}
