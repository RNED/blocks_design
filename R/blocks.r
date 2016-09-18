#' @title Block designs 
#' 
#' @description
#' 
#' Constructs randomized nested and crossed block designs for factorial and unstructured treatment sets where treatments can have any arbitrary levels of replication
#' and blocks can have any arbitrary feasible depth of nesting.
#' 
#' @details
#' 
#' The \code{blocksdesign} package constructs arbitrary block designs with arbitrary depth of nesting and arbitrary crossed row-and-column blocks
#' in each nested stratum. Where one classification of a row-and-column design has a single block level, the design automatically defaults to a simple nested blocks design
#' (if both classifications of a stratum have a single level that stratum is dropped from the design).       
#' 
#' The treatments design can be either an unstructured treatment set with an arbitrary number of treatments each with an arbitrary replication or a factorial
#' treatment design with an arbitrary combination of quantitative or qualitative level treatment factors.
#' 
#' For unstructured treatment sets, the \code{treatments} parameter must be a numeric vector partitioning the total number of treatments
#' into sets of equally replicated treatments. In this mode, the \code{replicates} parameter is a numeric vector giving
#' the required replication for each treatment set in the \code{treatments} partition. The \code{treatments} and \code{replicates} parameter vectors must be of equal length and the
#' sum of the cross products of the two parameter vectors must be the total number of plots.  
#' 
#' For factorial treatment sets, the \code{treatments} parameter must be a dataframe with a column of treatment levels for each factor of the design where the factorial design
#' comprises all possible combinations of the levels of the treatment factors. In this mode, the \code{replicates} parameter is a single number for the number of replications
#' of the full factorial design.
#' 
#' The default model for a fully crossed factorial \code{treatments} design is the full set of factorial effects but the \code{models} parameter can be used to define a model for 
#' any feasible reduced set of factorial effects. The \code{models} parameter must be a valid model formulation as defined by the \code{\link[stats]{model.matrix}} package
#' and any valid formulation for that package can be used to define a reduced set of factorial effects for the \code{blocksdesign} package. The \code{models} parameter is
#' ignored for any unstructured treatment sets.      
#' 
#' The \code{rows} vector, if specified, defines the nested row blocks in each nested stratum taken in order from the highest to the lowest.
#' The first number, if any, is the number of rows in the blocks of the top stratum, the second, if any, is the number of rows in the nested blocks of
#'  the second stratum, the third, if any, is the number of rows in the nested blocks of the third stratum and so on for all the required strata in the design. 
#'  
#' The \code{columns} vector, if specified, defines the nested column blocks in each nested stratum taken in order from the highest to the lowest. 
#' The first number, if any, is the number of columns in the blocks of the top stratum, the second, if any, is the number of columns in the nested blocks of
#'  the second stratum, the third, if any, is the number of columns in the nested blocks of the third stratum and so on for all the required strata in the design. 
#'   
#' The \code{rows} and \code{columns} vectors, if defined, must be of equal length and if a simple nested blocks design is required in 
#' any particular stratum the row or the column blocks classifications in that stratum should be set to unity.
#' 
#' If both the \code{rows} vector and the \code{columns} vector are null, the default block design will be a single set of orthogonal
#' main blocks equal in number to the highest common factor of the replication numbers. If the \code{rows} vector is defined but the \code{columns} vector
#' is null, the design will comprise simple nested blocks in each stratum defined by the \code{rows} vector.
#' 
#' Block sizes are always as nearly equal as possible and will never differ by more than a single plot in any particular classification.
#' Row blocks and column blocks must always contain at least two plots per block and this restriction will constrain the 
#' permitted numbers of rows and columns in the various strata of a design.
#' 
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number of single unreplicated treatments to give augmented blocks
#' that never differ in size by more than a single plot. General crossed block designs are more complex and currently 
#' the algorithm will only accommodate single unreplicated treatments in a crossed block design if the block sizes of the replicated part 
#' of the design are all equal in each stratum of the design.
#' 
#'  For row-and-column designs, the algorithm will unconditionally optimize the rows stratum and will then optimise the columns stratum conditional on the row blocks remaining 
#'  unchanged. Row-and-column designs contain nested blocks in the row-by-column intersections and these blocks may contain two or more plots. If the row-by-column intersection 
#'  blocks contain useful treatment information and the \code{weighted} option is TRUE (default) the algorithm will optimize the column design conditional on improving swaps 
#'  between columns also being improving swaps between blocks. If the \code{weighted} option is FALSE the column blocks are optimized ignoring blocks.           
#'  
#'  For 2 x 2 row-and-column designs with complete replicate rows and complete replicate columns, however, one treatment contrast will always be confounded
#'   with the row-by-column interaction and for these designs, it is impossible to nest blocks in the row-by-column intersections. 
#'  Instead, we recommend a simple nested blocks design with two complete or four incomplete main blocks. 
#'  
#'  Lattice designs where v is a prime-power require the \code{\link[crossdes]{MOLS}} package.
#' 
#'  The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order. \cr
#'  \item  A table showing the replication number of each treatment in the design. \cr
#'  \item  A table showing the block levels and the achieved D- and A-efficiency factors for each blocks stratum together with A-efficiency upper bounds, where available. \cr
#'  \item  Plans showing the allocation of treatments to blocks or to rows and to columns in the bottom stratum of the design. \cr
#' } 
#' 
#' @param treatments  either a dataframe with columns for the individual factors of the design or a vector giving a partition of the total required number of treatments into sets of equally replicated treatments.
#' 
#' @param replicates  either a single replication number for a factorial treatment set or a vector giving the replication number of each equally replicated treatment set in the \code{treatments} vector.
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
#' @param weighted  an option for row-and-column designs with two or more plot in the row.column blocks. \code{weighted}=TRUE
#' improves the row.columns blocks efficiency but may decrease the columns blocks efficiency. \code{weighted}=FALSE gives the best column
#' blocks efficiency but may give poor row.columns blocks efficiency. The default is weighted=TRUE.   
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
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. http://CRAN.R-project.org/package=crossdes
#' 
#' @examples
#' 
#' # Second-order model for two qualitative and two quntitative factors in 4 randomized blocks
#' TF=data.frame(F1=gl(2,36), F2=gl(3,12,72), V1=rep(rep(1:3,each=4),6), V2=rep(1:4,18))
#' model=" ~ F1*F2 + V1*V2 + I(V1^2) + I(V2^2) + F1:V1 + F1:V2 + F2:V1 + F2:V2"
#' blocks(treatments=TF,replicates=1,model=model,rows=4)
#' 
#' # 3 treatments x 2 replicates, 2 treatments x 4 replicates and 4 treatments x 3 replicates  
#' # the hcf of the replication numbers is 1 therefore the default design is completely randomized 
#' blocks(treatments=c(3,2,4),replicates=c(2,4,3))
#' 
#' # 4 treatments x 4 replicates with 2 main rows each containing two complete replicates  
#' blocks(treatments=4,replicates=4,rows=2)
#' 
#' # 50 treatments x 4 replicates with 4 main rows and 5 nested sub-rows in each main block 
#' blocks(treatments=50,replicates=4,rows=c(4,5))
#' 
#' # as above but with 20 single replicate treatments giving one extra treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),rows=c(4,5))
#' 
#' # 64 treatments x 2 replicates with 2 main rows and four succesively nested 2-level factors
#' blocks(treatments=64,replicates=2,rows=c(2,2,2,2,2))
#' 
#' # 64 treatments x 2 replicates with 2 main rows and five succesively nested 2-level factors
#' \dontrun{ blocks(treatments=64,replicates=2,rows=c(2,2,2,2,2,2)) }
#' 
#' # 6 replicates of 6 treatments in 4 rows of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,rows=4)
#' 
#' # 4 replicates of 13 treatments arranged in a 13 x 4 Youden rectangle 
#' blocks(treatments=13,replicates=4,rows=13,columns=4)
#' 
#' # 64 treatments x 2 replicates with nested 8 x 8 row-and-column designs in two main blocks 
#' blocks(treatments=64,replicates=2,rows=c(2,8),columns=c(1,8)) 
#' 
#' # 64 treatments x 2 replicates with two main blocks and a 4 x 4 row-and-column design  
#' # nested in each main block and weighted = TRUE by default.
#' blocks(treatments=64,replicates=2,rows=c(2,4),columns=c(1,4),searches=12) 
#' 
#' # 64 treatments x 2 replicates with two main blocks and a 4 x 4 row-and-column design  
#' # nested in each main block. weighted = FALSE
#' blocks(treatments=64,replicates=2,rows=c(2,4),columns=c(1,4),weighted=FALSE,searches=12)
#' 
#' # 2**9 treatments x 2 replicates in 2**9 blocks giving a fully saturated block design 
#' # (requires a considerable time to run!)
#' \dontrun{ d=blocks(2**9,2,rep(2,9),searches=1) }
#'          
#' @export
#' @importFrom stats anova lm
#' 
blocks = function(treatments,replicates,rows=HCF(replicates),columns=NULL,model=NULL,searches=NULL,seed=sample(10000,1),jumps=1,weighted=TRUE) { 
  options(contrasts=c('contr.sum','contr.poly'))

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
  # Tests a given number for primality and returns TRUE or FALSE
  # ********************************************************************************************************************************************************
  isPrime=function(v) {
    if (v < 4) return(TRUE) 
    if ( isTRUE(all.equal(v %% 2,0)) ||  isTRUE(all.equal(v %% 3,0)) ) return(FALSE) 
    if (v<25) return(TRUE)
    for(i in  6*seq_len(length(floor((sqrt(v)+1)/6)))        )
      if ( isTRUE(all.equal(v %% (i-1) , 0)) ||   isTRUE(all.equal(v %% (i+1) , 0)) ) return(FALSE) 
    return(TRUE)
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
  # Random swaps
  # ********************************************************************************************************************************************************    
  Swaps=function(TF,MF,BF,pivot,rank,n,Restrict) {
    candidates=NULL
    while (isTRUE(all.equal(length(candidates),0))) {
      if (rank<(n-1)) s1=sample(pivot[(1+rank):n],1) else s1=pivot[n]
      candidates = seq_len(n)[MF==MF[s1] & Restrict==Restrict[s1] & BF!=BF[s1] & TF!=TF[s1]]
    }
    if ( length(candidates)>1 )
      s2=sample(candidates,1) else s2=candidates[1] 
      s=c(s1,s2)
  }
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  NonSingular=function(TF,MF,BF,Restrict) { 
    BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
    BM[cbind(1:length(BF),BF)]=1
    fullrank=nlevels(TF)+nlevels(BF)-1
    TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
    TM[cbind(seq_len(length(TF)),TF)]=1
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    times=0
    n=length(TF)
    while (rank<fullrank & times<1000) {
      times=times+1
      s=Swaps(TF,MF,BF,pivot,rank,n,Restrict)
      rindex=seq_len(length(TF))
      rindex[c(s[1],s[2])]=rindex[c(s[2],s[1])]
      newQ=qr(t(cbind(BM,TM[rindex,])))
      if (isTRUE(all.equal(newQ$rank,rank)) || newQ$rank>rank) { 
        TF=TF[rindex]
        TM=TM[rindex,]
        rank=newQ$rank
        pivot=newQ$pivot
      } 
    }
    if (times>999) TF=NULL
    return(TF)
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
  DMax=function(MTT,MBB,MTB,Mtt,Mbb,Mtb,TF,weighted,Restrict,BF,Blocks){  
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
        if (!is.null(Mbb)&&isTRUE(weighted)) {
          Tb=Mtb[TF[s],Blocks[s],drop=FALSE]-tcrossprod(Mtb[cbind(TF[s],Blocks[s])],rep(1,nSamp[k]))
          dmat=(Tb+t(Tb)+1)**2-
            (2*Mtt[TF[s]  ,TF[s],  drop=FALSE]-tcrossprod(Mtt[cbind(TF[s]  ,TF[s])]  + rep(1,nSamp[k]) ) + tcrossprod(Mtt[cbind(TF[s],  TF[s])])   + 1)*
            (2*Mbb[Blocks[s],Blocks[s],drop=FALSE]-tcrossprod(Mbb[cbind(Blocks[s],Blocks[s])]+ rep(1,nSamp[k]) ) + tcrossprod(Mbb[cbind(Blocks[s],Blocks[s])]) + 1)
          is.na(dmat[ dmat<1|is.na(dmat) ] ) = TRUE
          dMat=dMat*dmat
        } 
        sampn=which.max(dMat)
        i=1+(sampn-1)%%nSamp[k]
        j=1+(sampn-1)%/%nSamp[k]
        maxd=dMat[i,j]
        if  (  maxd>1  && dMat[i,j]>1  &&  !isTRUE(all.equal(dMat[i,j],1))  ) {
          improved=TRUE
          relD=relD*maxd
          up=UpDate(MTT,MBB,MTB,TF[s[i]],TF[s[j]],BF[s[i]],BF[s[j]])
          MTT=up$MTT
          MBB=up$MBB
          MTB=up$MTB
          if (!is.null(Mbb)&&isTRUE(weighted)) {
            up=UpDate(Mtt,Mbb,Mtb,TF[s[i]],TF[s[j]],Blocks[s[i]],Blocks[s[j]])
            Mtt=up$MTT
            Mbb=up$MBB
            Mtb=up$MTB
          }
          TF[c(s[i],s[j])]=TF[c(s[j],s[i])]
        }
      } 
      if (improved) next
      if (sum(nSamp) < min(length(TF),512)) nSamp=pmin(mainSizes,2*nSamp) else break
    }
    list(MTT=MTT,MBB=MBB,MTB=MTB,Mtt=Mtt,Mbb=Mbb,Mtb=Mtb,TF=TF,relD=relD)
  }  
  
  # ******************************************************************************************************************************************************** 
  #  Searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # ********************************************************************************************************************************************************
  Optimise=function(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,Mtt,Mbb,Mtb,weighted)  {
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
    if (regReps && identical(max(breps),min(breps))  )
    bound=upper_bounds( length(TF), nlevels(TF), nlevels(BF) ) else bound=NA
    for (r in 1:searches) {
      dmax =  DMax(MTT,MBB,MTB,Mtt,Mbb,Mtb,TF,weighted,Restrict,BF,Blocks)
      if ( !isTRUE(all.equal(dmax$relD,1)) && dmax$relD>1) {
        relD=relD*dmax$relD
        TF=dmax$TF
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB 
        Mtt=dmax$Mtt
        Mbb=dmax$Mbb
        Mtb=dmax$Mtb 
        if (!isTRUE(all.equal(relD,globrelD)) && relD>globrelD) {
          globTF=TF
          globrelD=relD
          
          if ( !is.na(bound) && isTRUE(all.equal(bound,optEffics(globTF,BF)[2]))) break
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
          if (!is.null(Mbb)&&isTRUE(weighted)) 
            dswap = (1+Mtb[TF[s[1]],Blocks[s[2]]]+Mtb[TF[s[2]],Blocks[s[1]]]-Mtb[TF[s[1]],Blocks[s[1]]]-Mtb[TF[s[2]],Blocks[s[2]]])**2-
            (2*Mtt[TF[s[1]],TF[s[2]]]-Mtt[TF[s[1]],TF[s[1]]]-Mtt[TF[s[2]],TF[s[2]]])*(2*Mbb[Blocks[s[1]],Blocks[s[2]]]-Mbb[Blocks[s[1]],Blocks[s[1]]]-Mbb[Blocks[s[2]],Blocks[s[2]]])  
          else
            dswap=1
          if (Dswap>.1 & dswap>.1 | counter>1000) break
        }
        if (counter>1000) return(globTF) # no non-singular swaps
        relD=relD*Dswap 
        up=UpDate(MTT,MBB,MTB,TF[s[1]],TF[s[2]], BF[s[1]], BF[s[2]])
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        if (!is.null(Mbb)&&isTRUE(weighted)) {
          up=UpDate(Mtt,Mbb,Mtb,TF[s[1]],TF[s[2]],Blocks[s[1]],Blocks[s[2]])
          Mtt=up$MTT
          Mbb=up$MBB
          Mtb=up$MTB
        }
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
    if (is.factor(TF)) TF=NonSingular(TF,Main,Columns,Rows)
    if (is.null(TF)) return(TF)
    main=nlevels(Main)
    ncols=nlevels(Columns)/main
    nrows=nlevels(Rows)/main
    nblocks=nrows*ncols
    CM=Contrasts(Main,Columns)[,rep(c(rep(TRUE,(ncols-1)),FALSE),main),drop=FALSE]
    if (is.factor(TF)) TM=Contrasts(Main,TF)[,-nlevels(TF),drop=FALSE] else TM=TF
    V=chol2inv(chol(crossprod(cbind(CM,TM))))
    indicv=seq(ncol(CM)+1, ncol(TM)+ncol(CM))
    MTT=rbind(          cbind( V[indicv,indicv,drop=FALSE], rep(0,ncol(TM))) ,  rep(0,(ncol(TM)+1)))
    MBB=matrix(0,nrow=nlevels(Columns),ncol=nlevels(Columns))
    MTB=matrix(0,nrow=(ncol(TM)+1),ncol=nlevels(Columns))
    MBB[seq_len(ncol(CM)),seq_len(ncol(CM))]=V[seq_len(ncol(CM)),seq_len(ncol(CM)),drop=FALSE]
    MTB[seq_len(ncol(TM)),seq_len(ncol(CM))]=V[indicv,seq_len(ncol(CM)),drop=FALSE]
    perm=c(rbind(matrix(seq_len(ncol(CM)),nrow=ncols-1,ncol=main),seq_len(main)+ncol(CM)))
    MTB=MTB[,perm]
    MBB=MBB[perm,perm] 
    Blocks=as.factor((as.numeric(Main)-1)*nblocks + ((as.numeric(Rows)-1)%%nrows)*ncols + (as.numeric(Columns)-1)%%ncols + 1)
    BM=Contrasts(Main,Blocks)[,rep( c(rep(TRUE,(nblocks-1)),FALSE),main),drop=FALSE]
    if ( (nunits-1)>=(main*nrows*ncols+ncol(TM)) && qr(t(cbind(BM,TM)))$rank==(ncol(BM)+ncol(TM)) && isTRUE(weighted))  {
      V=chol2inv(chol(crossprod(cbind(BM,TM))))
      indicv=seq_len(ncol(TM))+ncol(BM)
      Mtt=rbind(cbind(V[indicv,indicv,drop=FALSE],rep(0,ncol(TM))) ,rep(0,(ncol(TM)+1) ))
      Mbb=matrix(0,nrow=main*nblocks,ncol=main*nblocks)
      Mtb=matrix(0,nrow=(ncol(TM)+1),ncol=main*nblocks)
      Mbb[seq_len(ncol(BM)),seq_len(ncol(BM))]=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
      Mtb[seq_len(ncol(TM)),seq_len(ncol(BM))]=V[indicv,seq_len(ncol(BM)),drop=FALSE]
      perm=c(rbind(matrix(seq_len(ncol(BM)),nrow=nblocks-1,ncol=main),seq_len(main)+ncol(BM)))
      Mtb=Mtb[,perm]
      Mbb=Mbb[perm,perm] 
      TF=Optimise(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,Mtt,Mbb,Mtb,weighted)
    } else 
      TF=Optimise(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,NULL,NULL,NULL,weighted)
    return(TF)
  }  
  # *******************************************************************************************************************************************************
  # Optimize the nested Blocks assuming a possible set of Main block constraints Initial randomized starting design. 
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ******************************************************************************************************************************************************** 
  rowsOpt=function(TF,Main,Rows,weighted) { 
      TF=NonSingular(TF,Main,Rows,rep(1,length(TF)))
      if (is.null(TF)) return(TF)
      BM=Contrasts(Main,Rows)[, rep(c(rep(TRUE,((nlevels(Rows)/nlevels(Main))-1)),FALSE),nlevels(Main)),drop=FALSE]
      if (is.factor(TF)) TM=Contrasts(Main,TF)[,-nlevels(TF),drop=FALSE] else TM=TF
      V=chol2inv(chol(crossprod(cbind(BM,TM))))
      indicv=seq_len(ncol(TM))+ncol(BM)
      MTT=rbind(cbind(V[indicv,indicv,drop=FALSE],rep(0,(nlevels(TF)-1))),rep(0,nlevels(TF)))
      MBB=matrix(0,nrow=nlevels(Rows),ncol=nlevels(Rows))
      MTB=matrix(0,nrow=nlevels(TF),ncol=nlevels(Rows))
      MBB[seq_len(ncol(BM)),seq_len(ncol(BM))]=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
      MTB[seq_len(ncol(TM)),seq_len(ncol(BM))]=V[indicv,seq_len(ncol(BM)),drop=FALSE]
      perm=c(rbind(matrix(seq_len(ncol(BM)),nrow=nlevels(Rows)/nlevels(Main)-1,ncol=nlevels(Main)),seq_len(nlevels(Main))+ncol(BM)))
      MTB=MTB[,perm]
      MBB=MBB[perm,perm] 
      TF=Optimise(TF,Main,Rows,NULL,NULL,MTT,MBB,MTB,NULL,NULL,NULL,weighted)
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
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  optEffics=function(TF,BF) { 
    k=nlevels(BF)
    if (k==1) return(c(1,1))
    if (nlevels(TF)<=k) 
      e=eigen( (diag(nlevels(TF))-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(nlevels(TF)-1)] else    
        e=c(rep(1,(nlevels(TF)-k)), 
            eigen((diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])  
      return(round(c(mean(e)*prod(e/mean(e))^(1/length(e)),1/mean(1/e)),6))
  }
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for block designs 
  # ********************************************************************************************************************************************************     
  BlockEfficiencies=function(Design) {
    strata=ncol(Design)-1
    effics=matrix(NA,nrow=strata,ncol=2)
    for (i in seq_len(strata))    
      effics[i,]=optEffics(Design$Treatments,Design[,i])  
    bounds=rep(NA,strata)
    if (max(replicates)==min(replicates)) 
      for (i in seq_len(strata))  
        if ( nunits%%nlevels(Design[,i])==0 ) 
          bounds[i]=upper_bounds(nunits,nlevels(TF),nlevels(Design[,i]) )
    names=NULL
    for (i in 1:strata)
      names=c(names,paste("Blocks",i))
    blocks=NULL
    for (i in 1:strata) 
      blocks=c( blocks,nlevels(Design[,i]))
    efficiencies=data.frame(cbind(names,blocks,effics,bounds))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
    efficiencies
  }
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for row-and-column designs 
  # ********************************************************************************************************************************************************     
  RowColEfficiencies=function(Design,Blocks) {
    strata=ncol(Blocks)-1
    effics=matrix(NA,nrow=(3*strata),ncol=2)
    for (i in seq_len(strata)) {
      effics[3*(i-1)+1,]=optEffics(Design$Treatments,Design[,2*(i-1)+1])  
      effics[3*(i-1)+2,]=optEffics(Design$Treatments,Design[,2*(i-1)+2])  
      effics[3*(i-1)+3,]=optEffics(Design$Treatments,Blocks[,i+1])  
    }
    bounds=rep(NA,(3*strata))
    if (max(replicates)==min(replicates)) {
      for (i in seq_len(strata))  {
        if (nunits%%nlevels(Design[,2*(i-1)+1])==0)
          bounds[3*(i-1)+1]=upper_bounds(nunits,nlevels(TF),nlevels(Design[,2*(i-1)+1]))
        if (nunits%%nlevels(Design[,2*(i-1)+2])==0)
          bounds[3*(i-1)+2]=upper_bounds(nunits,nlevels(TF),nlevels(Design[,2*(i-1)+2]))
        if (nunits%%nlevels(Blocks[,i+1])==0)
          bounds[3*(i-1)+3]=upper_bounds(nunits,nlevels(TF),nlevels(Blocks[,i+1]))
      }
    }
    names=NULL
    for (i in 1:strata)
      names=c(names,paste("Rows",i),paste("Columns",i),paste("Rows x Columns",i))
    blocks=NULL
    for (i in 1:strata) 
      blocks=c( blocks,nlevels(Design[,2*(i-1)+1]),nlevels(Design[,2*(i-1)+2]),nlevels(Blocks[,i+1]))
    efficiencies=data.frame(cbind(names,blocks,effics,bounds))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
  }
  

  # ******************************************************************************************************************************************************** 
  # Some validation checks
  # ********************************************************************************************************************************************************     
  Validate=function(treatments,replicates,rows,columns,seed,jumps,searches) {
    if (missing(treatments) | missing(replicates)) stop(" Treatments or replicates not defined ")   
    if (is.null(treatments) | is.null(replicates)) stop(" Treatments or replicates list is empty ")   
    if ((!is.data.frame(treatments))&&(any(!is.finite(treatments)) | any(!is.finite(replicates)) | any(is.nan(treatments)) | any(is.nan(replicates)))) stop("All inputs must be positive integers")
    if ((!is.data.frame(treatments))&& (length(treatments)!=length(replicates))) stop("The number of treatments sets and the number of replication numbers must match")
    if ((!is.data.frame(treatments))&&(sum(treatments)==1)) stop("Designs with only one treatment are not useful for comparative experiments")
    if ((is.data.frame(treatments))&&(nrow(treatments)==1)) stop("Designs with only one treatment are not useful for comparative experiments")
    if (anyNA(rows) ) stop(" NA rows values not allowed") 
    if (any(!is.finite(rows)) | any(is.nan(rows))) stop(" rows can contain only finite integers ")
    if (min(rows)<1) stop(" All row block parameters must be at least one ")
    if (anyNA(columns) ) stop(" NA columns values not allowed") 
    if (any(!is.finite(columns)) | any(is.nan(columns)) ) stop("All column block parameters must be finite integers ")
    if (length(columns)>0 && length(columns)!=length(rows)  ) stop("The rows parameter and the columns parameter must be the same length (same number of crossed-blocks strata)")
    if (length(columns)>0 && min(columns)<1) stop("All column block parameters must be at least one ")
    if (!is.null(searches)) {
      if (anyNA(searches) ) stop(" NA searches values not allowed") 
      if ( any(!is.finite(searches)) | any(is.nan(searches))) stop(" Searches must be a finite integer ") 
      if (searches<1) stop(" Repeats must be at least one ") 
    }  
    if (!is.null(jumps)) {
      if (anyNA(jumps) ) stop(" NA jumps values not allowed") 
      if ( !all(is.finite(jumps)) | !all(!is.nan(jumps))) stop(" jumps must be a finite integer ") 
      if (jumps<1)  stop(" Random jumps must be at least one ")   
      if (jumps>10)  stop(" Too many random jumps ") 
    }    
    if (!is.null(seed)) {
      if (anyNA(seed) ) stop(" NA seed values not allowed") 
      if (any(!is.finite(seed)) | any(is.nan(seed))) stop(" Seed must be a finite integer ") 
    } 
  } 
  
  # *******************************************************************************************************************************************************
  # Returns v-1 cyclically generated v x v Latin squares plus the rows array (first) and the columns array (last). If v is prime, the squares are MOLS  
  # ********************************************************************************************************************************************************  
  cMOLS=function(v) {
    mols=lapply(0:(v-1),function(z){do.call(rbind, lapply(0:(v-1), function(j){ (rep(0:(v-1))*z +j)%%v} ))})
    mols[[v+1]]=t(mols[[1]])
    mols=mols[c(1,length(mols),2:(length(mols)-1))]
    mols=lapply(1:length(mols),function(z){ (z-1)*v+mols[[z]] })
    mols
  }
  # *******************************************************************************************************************************************************
  # Tests for and constructs  balanced lattice designs
  # ******************************************************************************************************************************************************** 
  lattice=function(v,r) { 
    TF=rep(NA,v*v*r) 
    if ( r<4 || (isPrime(v) && r<(v+2)) ) {
      TF=rep(seq_len(v*v),r)[order(unlist(cMOLS(v))[1:(r*v*v)])]
    } else if (r<(v+2)  && (v*v)%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      TF[1:(2*v*v)]=c(seq_len(v*v),seq_len(v*v)[order(rep(0:(v-1),v))]   )
      if (r>2) {
        index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==(v*v))
        mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])	
        for (i in 1:(r-2)  )
          TF[ c( (v*v*(i-1)+1) : (v*v*i)) + 2*v*v  ] = seq_len(v*v)[order(as.numeric(mols[,,i]))]
      }
    } else if (v==10  && r==4) {
      square1=c(1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
                5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8)
      square2=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
                6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6)
      TF=c(seq_len(100),seq_len(100)[order(rep(0:9,10))],seq_len(100)[order(square1)],seq_len(100)[order(square2)])
    }
    TF=as.factor(TF)
  }
  # *******************************************************************************************************************************************************
  # Tests for balanced trojan designs and constructs available designs
  # ******************************************************************************************************************************************************** 
  trojan=function(r,k) { 
    TF=rep(NA,(r*r*k)) 
    if (isPrime(r)) { 
      for (z in 1:k)
        for (y in 0:(r-1)) 
          for (x in 0:(r-1)) 
            TF[(x + y*r)*k + z]=(y+x*z)%%r + (z-1)*r +1
    } else if ((r*r)%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==(r*r))
      mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])	
      for (i in 1:r)
        for (j in 1:r)
          for (x in 1:k)
            TF[x+(j-1)*k+(i-1)*k*r]=mols[i,j,x]+(x-1)*r
    }
    TF=as.factor(TF)
  }
  # ******************************************************************************************************************************************************** 
  # Main body of rows design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************     
  Validate(treatments,replicates,rows,columns,seed,jumps,searches) 
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
  if (!is.data.frame(treatments)) plots=sum(treatments[replicates>1]*replicates[replicates>1]) else plots=nrow(treatments)
  if (is.null(searches)) searches=1+10000%/%plots
  if ((!is.data.frame(treatments))&&(cumrows[strata]*2>plots & max(replicates)>1)) stop("Too many row blocks for the available plots  - every row block must contain at least two (replicated) treatments")
  if ((!is.data.frame(treatments))&&(cumcols[strata]*2>plots & max(replicates)>1)) stop("Too many column blocks for the available plots  - every column block must contain at least two (replicated) plots")
  if ((!is.data.frame(treatments))&&(cumblocks[strata+1]>plots & cumrows[strata]>1 & cumcols[strata]>1 )) stop("Too many blocks - every row-by-column intersection must contain at least one replicated plot")
  if ((!is.data.frame(treatments))&&( isTRUE( sum(treatments) < 2 ) )) stop(paste("The number of treatments must be at least two "))  
 
   if (!is.data.frame(treatments)&&(sum(treatments*replicates) < (cumrows[strata] + cumcols[strata] + sum(treatments)-2) ))
    stop(paste("The total number of plots is",sum(treatments*replicates), "whereas the total required number of model parameters is", 
               cumrows[strata] + cumcols[strata] + sum(treatments)-2)) 
  if ((!is.data.frame(treatments))&&(max(replicates)==2 && length(rows)>1 ))
    for (i in seq_len(length(rows)-1)) 
      if (rows[i]==2 && columns[i]==2) stop( paste("Cannot have nested sub-blocks within a 2 x 2 semi-Latin square - try a nested sub-blocks design within 4 main blocks"))
  set.seed(seed)
  if( is.data.frame(treatments)) {
    TM=model.matrix(  as.formula(model)  ,TF)
    TM=scale( TM[,c(2:ncol(TM))],center = TRUE, scale = FALSE)
    TM=TM[rep(seq_len(nrow(TM)), replicates), ]
    rownames(TM) = seq(length=nrow(TM)) 
    TM=as.matrix(TM)
    nunits=nrow(TM)
  }  
  if( !is.data.frame(treatments)) {  
  if (max(replicates)==1) {
    treatments=as.factor(sample(nlevels(TF)))
    Efficiencies=data.frame("Blocks 1","1",1,1,1)
    colnames(Efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    Design=data.frame(rep(1,nlevels(TF)), treatments)  
    Design[]=lapply(Design, as.factor)
    colnames(Design)=c("Blocks","Treatments")
    Treatments=data.frame(table(Design[,"Treatments"]))
    Treatments[]=lapply(Treatments, as.factor) 
    colnames(Treatments)=c("Treatments","Replicates")
    Plan=as.data.frame(cbind(1,"",matrix(treatments,nrow=1,ncol=length(treatments))))
    colnames(Plan)=c("Blocks 1","Plots",rep(1:nlevels(TF)))
    Plan[]  = lapply(Plan, as.factor)
    } else {
      # remove any single replicate treatments
    fulltreatments=treatments
    fullreplicates=replicates
    treatments=treatments[replicates>1]
    replicates=replicates[replicates>1]
    nunits=sum(treatments*replicates)
    rowcol=as.numeric(rbind(rows,columns))
    cpblocks=c(1,cumprod(rows*columns))
    cprowcol=cumprod(rowcol)
    isrowcol=max(columns)>1
    fDesign=do.call(cbind,lapply(1:(2*strata),function(i) {  gl(  rowcol[i],   cprowcol[2*strata]/cprowcol[i], cprowcol[2*strata]  )    })) -1 
    fBlocksInStrata=do.call(cbind,lapply(1:(strata+1),function(i) {  gl(  cpblocks[i], cpblocks[strata+1]/cpblocks[i], cpblocks[strata+1]  )    })) -1
    if (isrowcol)  fDesign=data.frame(do.call(cbind,lapply(1:(2*strata), function(r){fDesign[,r]+fBlocksInStrata[,(r-1)%/%2+1]*rowcol[r] }))+1)
    if (!isrowcol) fDesign=data.frame(do.call(cbind,lapply(1:strata,     function(r){fDesign[,2*r-1]+fBlocksInStrata[,r]*rowcol[2*r-1] }))+1)
    fDesign[]=lapply(fDesign, as.factor) 
    fBlocksInStrata=data.frame(fBlocksInStrata+1)
    fBlocksInStrata[]=lapply(fBlocksInStrata, as.factor) 
    #permBlocks will randomise whole blocks in nested strata or will randomize rows and columns in each nested stratum of a crossed design 
    tempDesign=data.frame( do.call(cbind,lapply(1:ncol(fDesign), function(r){ sample(nlevels(fDesign[,r]))[fDesign[,r]] })) , seq_len(nrow(fDesign)   ))
    permBlocks=tempDesign[ do.call(order, tempDesign), ][,ncol(tempDesign)]
    # design
    regBlocks=TRUE
    blocksizes=nunits
    for (i in seq_len(strata)) {
      blocksizes=Sizes(blocksizes,i) 
      if (max(blocksizes) > min(blocksizes)) regBlocks=FALSE
    }
    if ( min(fullreplicates)==1 && max(fullreplicates)>1 &&  max(columns)>1 && regBlocks==FALSE )  
      stop("The algorithm cannot deal with irregular row-and-column designs containing single replicate treatments ")
    Design  = data.frame(fDesign[rep(seq_len(nrow(fDesign)),  blocksizes ),])  
    BlocksInStrata  = data.frame(fBlocksInStrata[rep(seq_len(nrow(fBlocksInStrata)),  blocksizes ),]) 
    Design[]= lapply(Design, as.factor) 
    BlocksInStrata[]= lapply(BlocksInStrata, as.factor) 
    hcf=HCF(replicates)
    sqrLattice=FALSE 
    fulltrojan=FALSE
    regReps=identical(length(replicates),as.integer(1))
    orthoMain=(regReps && (replicates[1]==rows[1]))
    v=sqrt(sum(treatments))  # dimension of a lattice square
    k=nunits/prod(rows*columns)  # block size 
    r=replicates[1]
    if (regReps && regBlocks && orthoMain && !isrowcol && identical(v,floor(v)) && identical(k,v) && identical(length(rows),as.integer(2))) {
    TF=lattice(v,r)
    if (!all(is.na(TF))) sqrLattice=TRUE
    }
    # given s orthogonal Latin squares of dimension r x r there are r x kr Trojan designs for r replicates of kr treatments in blocks of size k where k<=s
    if (regReps && regBlocks && orthoMain && isrowcol && identical(columns[1],r) && identical(length(rows),as.integer(1)) && identical(length(columns),as.integer(1)) && (k<r)) {
    TF=trojan(r,k)
    if (!all(is.na(TF))) fulltrojan=TRUE
    }
    # Treatment factors and levels ignoring any single replicate treatments
    if ((!sqrLattice) && (!fulltrojan)) { 
    trtReps=rep(fullreplicates,fulltreatments)
    TrtLabels=rep(1:sum(fulltreatments))[trtReps>1]
    for ( z in seq_len(10)) {
      TF=rep(rep(TrtLabels,trtReps[trtReps>1]/hcf),hcf)
      rand=sample(nunits)      
      TF=as.factor( TF[rand][order(rep(seq_len(hcf),each=nunits/hcf)[rand])] )
      for ( i in seq_len(strata)) {
        if (!isrowcol && rows[i]>1 && !all(hcf%%cumprod(rows[1:i])==0)) TF=rowsOpt(TF,BlocksInStrata[,i],Design[,i],weighted)
        if (is.null(TF)) break
        if (isrowcol &&  rows[i]>1 && !all(hcf%%cumprod(rows[1:i])==0)) TF=rowsOpt(TF,BlocksInStrata[,i],Design[,2*i-1],weighted)
        if (is.null(TF)) break
        if (isrowcol && columns[i]>1) TF=colsOpt(TF,BlocksInStrata[,i],Design[,2*i-1],Design[,2*i],weighted)
        if (is.null(TF)) break
      }
      if (!is.null(TF)) break
    }
    if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
    }
    # add back single rep treatments for nested stratum blocks only
    if ( min(fullreplicates)==1 && max(fullreplicates)>1  && ( max(columns)==1 || regBlocks==TRUE ) ) {
      replicates=fullreplicates
      treatments=fulltreatments
      fullreps=rep(replicates,treatments)
      nunits=sum(treatments*replicates)
      fullblocksizes=nunits
      for (i in seq_len(strata))
        fullblocksizes=Sizes(fullblocksizes,i)
      sblocksizes=fullblocksizes-blocksizes
      singleTF=     split(rep(1:sum(treatments))[fullreps==1],rep( 1:length(sblocksizes),sblocksizes))
      TrtsInBlocks= split(as.numeric(levels(TF))[TF],         rep(1:length(blocksizes),blocksizes))
      addBlocks=which(sblocksizes>0)
      for (i in 1:length(addBlocks)) 
        TrtsInBlocks[[addBlocks[i]]]=sample(append( TrtsInBlocks[[addBlocks[i]]] ,  singleTF[[i]]))
      TF=as.factor(unlist(TrtsInBlocks))
      blocksizes=fullblocksizes
      Design  = data.frame(fDesign[rep(seq_len(nrow(fDesign)),  blocksizes ),])  
      BlocksInStrata  = data.frame(fBlocksInStrata[rep(seq_len(nrow(fBlocksInStrata)),  blocksizes ),]) 
      Design[]= lapply(Design, as.factor) 
      BlocksInStrata[]= lapply(BlocksInStrata, as.factor) 
    }
    # Randomize
    D=as.data.frame(cbind(rep(1:length(blocksizes),blocksizes),sample(seq_len(nunits)),TF))
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
    # efficiencies
    if (isrowcol)     
      Efficiencies=RowColEfficiencies(Design,BlocksInStrata)
    else
      Efficiencies=BlockEfficiencies(Design)
    row.names(Efficiencies)=NULL
    # omit single level row or column strata in row and column designs
    if (isrowcol) Design[c(which(as.numeric(rbind(rows,columns))==1))]= list(NULL) 
    if (isrowcol) Plan[c(which(as.numeric(rbind(rows,columns))==1 ))]= list(NULL) 
    }
  }
  
  # treatment replications
  TreatmentsTable=data.frame(table(Design[,"Treatments"]))
  TreatmentsTable[]=lapply(TreatmentsTable, as.factor) 
  colnames(TreatmentsTable)=c("Treatments","Replicates")
  
  list(Treatments=TreatmentsTable,Efficiencies=Efficiencies,Plan=Plan,Design=Design,Seed=seed,Searches=searches,Jumps=jumps) 
} 