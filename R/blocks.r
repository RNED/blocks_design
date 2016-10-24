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
#' of the full factorial design. Only designs with a single blocks stratum, possibly resolved into complete replicate blocks for replicated treatment sets, 
#' are available for factorial designs.   
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
#' First-order model for five qualitative two level factors in 4 randomized blocks
#' TF=data.frame( F1=gl(2,16), F2=gl(2,8,32),  F3=gl(2,4,32), F4=gl(2,2,32) , F5=gl(2,1,32)  )
#' model=" ~ (F1+F2+F3+F4+F5)"
#' blocks(treatments=TF,replicates=1,model=model,rows=4)
#' 
#' # Second-order model for five qualitative two level factors in 4 randomized blocks
#' TF=data.frame( F1=gl(2,16), F2=gl(2,8,32),  F3=gl(2,4,32), F4=gl(2,2,32) , F5=gl(2,1,32)  )
#' model=" ~ (F1+F2+F3+F4+F5)*(F1+F2+F3+F4+F5)"
#' blocks(treatments=TF,replicates=1,model=model,rows=4)
#' blocks(treatments=TF,replicates=1,model=model,rows=c(2,2))
#' 
#' # Second-order design for four qualitative three level factors in 9 randomized blocks
#' TF=data.frame( F1=gl(3,27), F2=gl(3,9,81),  F3=gl(3,3,81), F4=gl(3,1,81)  )
#' model=" ~ (F1+F2+F3+F4)*(F1+F2+F3+F4)" # main effects and 2-factor interactions
#' blocks(treatments=TF,replicates=1,model=model,rows=9)
#' 
#' # Second-order model for two qualitative and two quantitative factors in 4 randomized blocks
#' TF=data.frame(F1=gl(2,36), F2=gl(3,12,72), V1=rep(rep(1:3,each=4),6), V2=rep(1:4,18))
#' model=" ~ F1*F2 + V1*V2 + I(V1^2) + I(V2^2) + F1:V1 + F1:V2 + F2:V1 + F2:V2"
#' blocks(treatments=TF,replicates=1,model=model,rows=4)
#' 
#' Explicit factorial model for 4 replicates of 13 treatments arranged in a 13 x 4 Youden rectangle 
#' TF=data.frame( Treatments=as.factor(1:13) )
#' model=" ~ Treatments"
#' blocks(treatments=TF,replicates=4,model=model,rows=13)
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
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
  # ********************************************************************************************************************************************************
  factUpDate=function(MTT,MBB,MTB,t,b) {
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
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
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
        s=sort(sample( seq_len(nrow(TF)) [Restrict==levels(Restrict)[k]], nSamp[k])) 
        TB=MTB[TF[s,1],BF[s],drop=FALSE]-tcrossprod(MTB[cbind(TF[s,1],BF[s])],rep(1,nSamp[k]))
        dMat=(TB+t(TB)+1)**2-
          (2*MTT[TF[s,1],TF[s,1],drop=FALSE]- tcrossprod(MTT[cbind(TF[s,1],TF[s,1])]+rep(1,nSamp[k])) + tcrossprod(MTT[cbind(TF[s,1],TF[s,1])]) + 1)*
          (2*MBB[BF[s],BF[s],drop=FALSE]- tcrossprod(MBB[cbind(BF[s],BF[s])]+rep(1,nSamp[k])) + tcrossprod(MBB[cbind(BF[s],BF[s])]) + 1)
        is.na(dMat[ dMat<0.00001|is.na(dMat) ] ) = TRUE
        if (!is.null(Mbb)&&isTRUE(weighted)) {
          Tb=Mtb[TF[s,1],Blocks[s],drop=FALSE]-tcrossprod(Mtb[cbind(TF[s,1],Blocks[s])],rep(1,nSamp[k]))
          dmat=(Tb+t(Tb)+1)**2-
            (2*Mtt[TF[s,1],TF[s,1],drop=FALSE]-tcrossprod(Mtt[cbind(TF[s,1],TF[s,1])] + rep(1,nSamp[k]) ) + tcrossprod(Mtt[cbind(TF[s,1], TF[s,1])]) + 1)*
            (2*Mbb[Blocks[s],Blocks[s],drop=FALSE]-tcrossprod(Mbb[cbind(Blocks[s],Blocks[s])]+ rep(1,nSamp[k]) ) + tcrossprod(Mbb[cbind(Blocks[s],Blocks[s])]) + 1)
          is.na(dmat[ dmat<tol|is.na(dmat) ] ) = TRUE
          dMat=dMat*dmat
        } 
        sampn=which.max(dMat)
        i=1+(sampn-1)%%nSamp[k]
        j=1+(sampn-1)%/%nSamp[k]
        if  (dMat[i,j]>tol) {
          improved=TRUE
          relD=relD*dMat[i,j]
          up=UpDate(MTT,MBB,MTB,TF[s[i],1],TF[s[j],1],BF[s[i]],BF[s[j]])
          MTT=up$MTT
          MBB=up$MBB
          MTB=up$MTB
          if (!is.null(Mbb)&&isTRUE(weighted)) {
            up=UpDate(Mtt,Mbb,Mtb,TF[s[i],1],TF[s[j],1],Blocks[s[i]],Blocks[s[j]])
            Mtt=up$MTT
            Mbb=up$MBB
            Mtb=up$MTB
          }
          TF[c(s[i],s[j]),1]=TF[c(s[j],s[i]),1]
        }
      } 
      
      if (improved) next
      if (sum(nSamp) < min(nrow(TF),512)) nSamp=pmin(mainSizes,2*nSamp) else break
    }
    list(MTT=MTT,MBB=MBB,MTB=MTB,Mtt=Mtt,Mbb=Mbb,Mtb=Mtb,TF=TF,relD=relD,TM=NULL)
  }  
  # ******************************************************************************************************************************************************** 
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
  # ********************************************************************************************************************************************************
  factDMax=function(MTT,MBB,MTB,Mtt,Mbb,Mtb,TF,weighted,Main,BF,Blocks,TM,BM,RCM) {  
    relD=1
    mainSizes=tabulate(Main)
    nSamp=pmin( rep(8,nlevels(Main)), mainSizes)
    repeat {
      improved=FALSE
      for (k in 1: nlevels(Main)) {
        s=sort(sample( seq_len(nrow(TF)) [Main==levels(Main)[k]], nSamp[k])) 
        TMB=crossprod(t(crossprod(t(TM[s,]),MTB)),t(BM[s,]))
        TMT=crossprod(t(crossprod(t(TM[s,]),MTT)),t(TM[s,]))
        BMB=crossprod(t(crossprod(t(BM[s,]),MBB)),t(BM[s,]))
        TMB=sweep(TMB,1,diag(TMB))
        TMT=sweep(TMT,1,diag(TMT))
        BMB=sweep(BMB,1,diag(BMB))
        dMat=(1+TMB+t(TMB))**2 - (TMT + t(TMT))*(BMB + t(BMB))
        is.na(dMat[ dMat<0.00001|is.na(dMat) ] ) = TRUE
        sampn=which.max(dMat)
        i=1+(sampn-1)%%nrow(dMat)
        j=1+(sampn-1)%/%nrow(dMat)
        if  ( dMat[i,j]>tol) {
          improved=TRUE
          relD=relD*dMat[i,j]
          t=TM[s[i],]-TM[s[j],]
          b=BM[s[j],]-BM[s[i],]
          up=factUpDate(MTT,MBB,MTB,t,b)
          MTT=up$MTT
          MBB=up$MBB
          MTB=up$MTB
          # unfinisfhed yet to be done
          if (!is.null(Mbb)&&isTRUE(weighted)) {
            t=TM[s[i],]-TM[s[j],]
            b=BM[s[j],]-BM[s[i],]
            up=factUpDate(MTT,MBB,MTB,t,b)
            up=UpDate(Mtt,Mbb,Mtb,TF[s[i]],TF[s[j]],Blocks[s[i]],Blocks[s[j]])
            Mtt=up$MTT
            Mbb=up$MBB
            Mtb=up$MTB
          }
          TF[c(s[i],s[j]),]=TF[c(s[j],s[i]),]
          TM[c(s[i],s[j]),]=TM[c(s[j],s[i]),]
        }
      } 
      if (improved) next
      if (sum(nSamp) < min(nrow(TF),512)) nSamp=pmin(mainSizes,2*nSamp) else break
    }
    list(MTT=MTT,MBB=MBB,MTB=MTB,Mtt=Mtt,Mbb=Mbb,Mtb=Mtb,TF=TF,relD=relD,TM=TM)
  }  
  # ******************************************************************************************************************************************************** 
  #  Searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # ********************************************************************************************************************************************************
  Optimise=function(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,Mtt,Mbb,Mtb,weighted,TM,BM,RCM)  {
    if (is.null(Rows)) {
      Restrict=rep(1,length(Main)) 
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
    if ( regReps && identical(max(breps),min(breps))  && !is.data.frame(treatments))
      bound=upper_bounds( nrow(TF), nlevels(TF[,1]), nlevels(BF) ) else bound=NA
    for (r in 1:searches) {
      if (is.data.frame(treatments)) 
        dmax =factDMax(MTT,MBB,MTB,Mtt,Mbb,Mtb,TF,weighted,Restrict,BF,Blocks,TM,BM,RCM) 
      else
        dmax =  DMax(MTT,MBB,MTB,Mtt,Mbb,Mtb,TF,weighted,Restrict,BF,Blocks)
      if (dmax$relD>tol) {
        relD=relD*dmax$relD
        TF=dmax$TF
        TM=dmax$TM
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB 
        Mtt=dmax$Mtt
        Mbb=dmax$Mbb
        Mtb=dmax$Mtb 
        if (!isTRUE(all.equal(relD,globrelD)) && relD>globrelD) {
          globTF=TF
          globrelD=relD
          if ( !is.na(bound) && isTRUE(all.equal(bound,optEffics(globTF[,1],BF)[2]))) break
        }
      }

      if (r==searches) break
      for (iswap in 1:jumps) {
        counter=0
        repeat {  
          counter=counter+1
          s1=sample(seq_len(nunits),1)
          altTF=apply(  sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,prod)==0 
          z= seq_len(nunits)[ Main==Main[s1] & Restrict==Restrict[s1] & BF!=BF[s1] & altTF==TRUE ]     
          if (length(z)==0) next
          if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z[1])

          if (is.data.frame(treatments)) {
            TMB=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],MTB)),BM[s[2],]-BM[s[1],] )
            TMT=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],MTT)),TM[s[2],]-TM[s[1],])
            BMB=crossprod(t(crossprod(BM[s[1],]-BM[s[2],],MBB)),BM[s[2],]-BM[s[1],])
            Dswap=(1+TMB)**2-TMT*BMB
          } else 
            Dswap = (1+MTB[TF[s[1],],BF[s[2]]]+MTB[TF[s[2],],BF[s[1]]]-MTB[TF[s[1],],BF[s[1]]]-MTB[TF[s[2],],BF[s[2]]])**2-
              (2*MTT[TF[s[1],],TF[s[2],]]-MTT[TF[s[1],],TF[s[1],]]-MTT[TF[s[2],],TF[s[2],]])*(2*MBB[BF[s[1]],BF[s[2]]]-MBB[BF[s[1]],BF[s[1]]]-MBB[BF[s[2]],BF[s[2]]])  

          if (!is.null(Mbb)&&isTRUE(weighted)) 
            dswap = (1+Mtb[TF[s[1],],Blocks[s[2]]]+Mtb[TF[s[2],],Blocks[s[1]]]-Mtb[TF[s[1],],Blocks[s[1]]]-Mtb[TF[s[2],],Blocks[s[2]]])**2-
            (2*Mtt[TF[s[1],],TF[s[2],]]-Mtt[TF[s[1],],TF[s[1],]]-Mtt[TF[s[2],],TF[s[2],]])*(2*Mbb[Blocks[s[1]],Blocks[s[2]]]-Mbb[Blocks[s[1]],Blocks[s[1]]]-Mbb[Blocks[s[2]],Blocks[s[2]]])  
          else
            dswap=1
          if (Dswap>.1 & dswap>.1 | counter>1000) break
        }

        if (counter>1000) return(globTF) # no non-singular swaps
        relD=relD*Dswap 
        if (is.data.frame(treatments)) {
          t=TM[s[1],]-TM[s[2],]
          b=BM[s[2],]-BM[s[1],]
          up=factUpDate(MTT,MBB,MTB,t,b)
        } else {
          up=UpDate(MTT,MBB,MTB,TF[s[1],],TF[s[2],], BF[s[1]], BF[s[2]])
        }
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        
        if (!is.null(Mbb)&&isTRUE(weighted)) {
          up=UpDate(Mtt,Mbb,Mtb,TF[s[1],],TF[s[2],],Blocks[s[1]],Blocks[s[2]])
          Mtt=up$MTT
          Mbb=up$MBB
          Mtb=up$MTB
        }
        TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]  
        TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]  
      } 
    }
    globTF
  } 
  
  # ******************************************************************************************************************************************************** 
  # Contrasts for factor NF centered within the levels of factor MF to ensure that NF information is estimated within the levels of factor MF only  
  # ********************************************************************************************************************************************************
  Contrasts=function(MF,NF) {
    NM=matrix(0,nrow=length(NF),ncol=nlevels(NF))
    NM[cbind(seq_len(length(NF)),NF)]=1 # factor indicator matrix  
    NM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(NM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
  }
  
  # *******************************************************************************************************************************************************
  # Optimize the nested Blocks assuming a possible set of Main block constraints Initial randomized starting design. 
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ******************************************************************************************************************************************************** 
  rowsOpt=function(TF,MF,BF) { 
    if (is.data.frame(treatments)) {  
    TM=model.matrix(as.formula(model),TF)[,-1,drop=FALSE] # drops mean contrast
    TM=do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(TM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
    } else TM=Contrasts(MF,TF[,1])[,-nlevels(TF[,1]),drop=FALSE] # drops last treatment contrast
    BM=Contrasts(MF,BF)[, -seq( nlevels(BF)/nlevels(MF), nlevels(BF) , by=nlevels(BF)/nlevels(MF) ) ,drop=FALSE]
    nonsing=NonSingular(TF,MF,BF,TM,BM)
    TF=nonsing$TF
    TM=nonsing$TM
    if (is.null(TF)) return(TF)
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    indicv=seq( from=(ncol(BM)+1), to=(ncol(BM)+ncol(TM))  )
    if (!is.data.frame(treatments)) {
      MTT=rbind(  cbind( V[indicv,indicv,drop=FALSE],rep(0,length(indicv))) , rep(0,(1+length(indicv))) )
      MBB=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
      MBB[seq_len(ncol(BM)),seq_len(ncol(BM))]=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
      MTB=matrix(0,nrow=nlevels(TF[,1]),ncol=nlevels(BF))
      MTB[seq_len(ncol(TM)),seq_len(ncol(BM))]=V[indicv,seq_len(ncol(BM)),drop=FALSE]
      reorder=c(rbind(matrix(seq_len(ncol(BM)),nrow=nlevels(BF)/nlevels(MF)-1,ncol=nlevels(MF)),seq(ncol(BM)+1, nlevels(BF))))
      MBB=MBB[reorder,reorder] 
      MTB=MTB[,reorder]
    } else {
      MBB=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
      MTB=V[indicv,seq_len(ncol(BM)),drop=FALSE]
      MTT=V[indicv,indicv,drop=FALSE] 
    }
    TF=Optimise(TF,MF,BF,NULL,NULL,MTT,MBB,MTB,NULL,NULL,NULL,weighted,TM,BM,NULL)
    if (is.null(TF))  stop(" Unable to find a solution design in the rows stratum ") 
    TF
  } 
  # ******************************************************************************************************************************************************** 
  # Optimize the nested columns blocks assuming a possible set of Main block constraints Initial randomized starting design. 
  # If the row-column intersections contain 2 or more plots a weighted (columns + w*rows.columns) model is fitted for w>=0 and w<1 
  # ********************************************************************************************************************************************************    
  colsOpt=function(TF,model,Main,Rows,Columns,Blocks,weighted,TM) { 
    if (!is.data.frame(treatments)) 
      TM=Contrasts(Main,TF) # treatments nested in overall mean
    TM=TM[,-ncol(TM),drop=FALSE] # drops last treatment contrast
    CM=Contrasts(Main,Columns)[, -c(1:nlevels(Main))*(nlevels(Columns)/nlevels(Main)) ,drop=FALSE]
    nonsing=NonSingular(TF,Rows,Columns,TM,CM)
    TF=nonsing$TF
    TM=nonsing$TM
    if (is.null(TF)) return(TF)
    V=chol2inv(chol(crossprod(cbind(CM,TM))))
    main=nlevels(Main)
    ncols=nlevels(Columns)/main
    nrows=nlevels(Rows)/main
    nblocks=nrows*ncols
    indicv=seq( (ncol(CM)+1), (ncol(CM)+ncol(TM))  )
    if (!is.data.frame(treatments)) {
      perm=c(rbind(matrix(seq_len(ncol(CM)), nrow=ncols-1, ncol=main), seq((ncol(CM)+1),(main+ncol(CM)))))
      MTT=matrix(0,nrow=nlevels(TF[,1]),ncol=nlevels(TF[,1]))
      MTT[seq_len(ncol(TM)),seq_len(ncol(TM))]=V[indicv,indicv,drop=FALSE]
      MBB=matrix(0,nrow=nlevels(Columns),ncol=nlevels(Columns))
      MBB[seq_len(ncol(CM)),seq_len(ncol(CM))]=V[seq_len(ncol(CM)),seq_len(ncol(CM)),drop=FALSE]
      MBB=MBB[perm,perm] 
      MTB=matrix(0,nrow=nlevels(TF[,1]),ncol=nlevels(Columns))
      MTB[seq_len(ncol(TM)),seq_len(ncol(CM))]=V[indicv,seq_len(ncol(CM)),drop=FALSE]
      MTB=MTB[,perm]
    } else {
      MBB=V[seq_len(ncol(CM)),seq_len(ncol(CM)),drop=FALSE]
      MTB=V[indicv,seq_len(ncol(CM)),drop=FALSE]
      MTT=V[indicv,indicv,drop=FALSE] 
    }  
    BM=Contrasts(Main,Blocks)
    T=table(Main,Blocks)
    BM=BM[,-unlist(lapply(1:nrow(T),function(i) {which(T[i,] > 0)[1]})),drop=FALSE] # drops first contrast in each nested block
    
    indicv=seq( (ncol(BM)+1), (ncol(BM)+ncol(TM))  )
    if ( nunits>=(main*nrows*ncols+ncol(TM)+1) && qr(t(cbind(BM,TM)))$rank==(ncol(BM)+ncol(TM)) && isTRUE(weighted))  {
      V = chol2inv(chol(crossprod(cbind(BM,TM))))
      if (!is.data.frame(treatments)) {
        perm=c(rbind(matrix(seq_len(ncol(BM)),nrow=nblocks-1,ncol=main),     seq((ncol(BM)+1),(main+ncol(BM))) ))
        Mtt = matrix(0,nrow=nlevels(TF[,1]),ncol=nlevels(TF[,1]))
        Mtt[seq_len(ncol(TM)),seq_len(ncol(TM))]=V[indicv,indicv,drop=FALSE]
        Mbb=matrix(0,nrow=nlevels(Blocks),ncol=nlevels(Blocks))
        Mbb[seq_len(ncol(BM)),seq_len(ncol(BM))]=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
        Mbb=Mbb[perm,perm] 
        Mtb=matrix(0,nrow=(ncol(TM)+1),ncol=main*nblocks)
        Mtb[seq_len(ncol(TM)),seq_len(ncol(BM))]=V[indicv,seq_len(ncol(BM)),drop=FALSE]
        Mtb=Mtb[,perm]
        Mbb=Mbb[perm,perm] 
      } else {
        Mbb=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
        Mtb=V[indicv,seq_len(ncol(BM)),drop=FALSE]
        Mtt=V[indicv,indicv,drop=FALSE] 
      }
      TF=Optimise(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,Mtt,Mbb,Mtb,weighted,TM,CM,BM)
    } else {
      TF=Optimise(TF,Main,Rows,Columns,Blocks,MTT,MBB,MTB,NULL,NULL,NULL,weighted,TM,CM,BM)
    }
    if (is.null(TF))  stop(" Unable to find a solution design in the columns stratum ") 
    TF
  }  
 
  # ******************************************************************************************************************************************************** 
  # Random swaps
  # ********************************************************************************************************************************************************    
  Swaps=function(TF,MF,BF,pivot,rank,nunits) {
    candidates=NULL
    while (isTRUE(all.equal(length(candidates),0))) {
      if (rank<(nunits-1)) s1=sample(pivot[(1+rank):nunits],1) else s1=pivot[nunits]
      altTF=apply(  sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,prod)==0
      candidates = seq_len(nunits)[ MF==MF[s1] & BF!=BF[s1] & altTF==TRUE ]
    }
    if ( length(candidates)>1 ) s2=sample(candidates,1) else s2=candidates[1] 
    return(c(s1,s2))
  }
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  NonSingular=function(TF,MF,BF,TM,BM) {
    fullrank=ncol(BM)+ncol(TM)
    
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    times=0
    while (rank<fullrank & times<1000) {
      times=times+1
      s=Swaps(TF,MF,BF,pivot,rank,nunits)
      tindex=seq_len(length(BF))
      tindex[c(s[1],s[2])]=tindex[c(s[2],s[1])]
      newQ=qr(t(cbind(BM,TM[tindex,])))
      if (newQ$rank>rank) {   
        TF=TF[tindex,,drop=FALSE]
        TM=TM[tindex,,drop=FALSE]
        rank=newQ$rank
        pivot=newQ$pivot
      } 
    }
    if (times>999) TF=NULL
    list(TF=TF,TM=TM)
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
  # Finds efficiency factors for block designs 
  # ********************************************************************************************************************************************************     
  BlockEfficiencies=function(Design) {
    effics=matrix(NA,nrow=strata,ncol=2)
    for (i in seq_len(strata))    
      effics[i,]=optEffics(Design[,ncol(Design)],Design[,i])  
    
    bounds=rep(NA,strata)
    if (max(replicates)==min(replicates)) 
      for (i in seq_len(strata))  
        if ( nunits%%nlevels(Design[,i])==0 ) 
          bounds[i]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,i]) )
  
    names =unlist(lapply(1:strata, function(j) {paste0("Stratum_",j)}))
    blocks=unlist(lapply(1:strata, function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocks,effics,bounds))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
  }
  
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for row-and-column designs 
  # ********************************************************************************************************************************************************     
  RowColEfficiencies=function(Design) {
    effics=matrix(NA,nrow=(3*strata),ncol=2)
    for (i in seq_len(strata)) {
      effics[3*(i-1)+1,]=optEffics(Design[,ncol(Design)],Design[,3*(i-1)+1])  
      effics[3*(i-1)+2,]=optEffics(Design[,ncol(Design)],Design[,3*(i-1)+2])  
      effics[3*(i-1)+3,]=optEffics(Design[,ncol(Design)],Design[,3*(i-1)+3])  
    }
    bounds=rep(NA,(3*strata))
    if (max(replicates)==min(replicates)) {
      for (i in seq_len(strata))  {
        if (nunits%%nlevels(Design[,3*(i-1)+1])==0)
          bounds[3*(i-1)+1]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+1]))
        if (nunits%%nlevels(Design[,3*(i-1)+2])==0)
          bounds[3*(i-1)+2]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+2]))
        if (nunits%%nlevels(Design[,3*(i-1)+3])==0)
          bounds[3*(i-1)+3]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+3]))
      }
    }
    names =unlist(lapply(1:strata, function(j) {c(paste("Rows",j),paste("Columns",j),paste("Rows x Columns",j))}))
    blocks=unlist(lapply(1:strata, function(j) {c( nlevels(Design[,3*(j-1)+1]),nlevels(Design[,3*(j-1)+2]),nlevels(Design[,3*(j-1)+3]))}))
    efficiencies=data.frame(cbind(names,blocks,effics,bounds))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
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
  # Calculates D and A-efficiency factors for treatment factors TF assuming block factor BF
  # ********************************************************************************************************************************************************
  factEffics=function(Design) { 
    TF=Design[, c( (ncol(Design)-ncol(TF)+1):ncol(Design)),drop=FALSE]
    TX=model.matrix(as.formula(model),TF)[,-1,drop=FALSE] # drops mean contrast
    names =unlist(lapply(1:strata, function(j) {paste0("Stratum_",j)}))
    blocks=unlist(lapply(1:strata, function(j) {nlevels(Design[,j])}))
    effics=NULL
    Design=data.frame(as.factor(rep(1,nrow(Design))),Design)
    for (i in seq_len(strata)) { 
      MF=Design[,i]
      BF=Design[,i+1]
      TM=do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(TX[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
      BM=Contrasts(MF,BF)[, -seq( nlevels(BF)/nlevels(MF), nlevels(BF) , by=nlevels(BF)/nlevels(MF) ) ,drop=FALSE]
      TB=crossprod(TM,BM)
      RI=backsolve(  chol(crossprod(TM)) ,diag(ncol(TM)))
      QI=backsolve(chol(crossprod(BM)),diag(ncol(BM)))
      U=crossprod(t(crossprod(RI,TB)),QI)
      effics=round(c(effics,det( diag(ncol(TM))-tcrossprod(U))**(1/ncol(TM))),4)
    }
    efficiencies=data.frame(cbind(names,blocks,effics))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies")
    efficiencies
  }
  # ******************************************************************************************************************************************************** 
  # Efficiency factors for unreplicated randomized designs 
  # ********************************************************************************************************************************************************     
  UnrepEfficiencies=function() {
    efficiencies=data.frame(cbind("Stratum_1",1,1,1,1))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
  }
  
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
  # Contrasts for factor NF centered within the levels of factor MF to ensure that NF information is estimated within the levels of factor MF only  
  # ********************************************************************************************************************************************************
  trtMat=function(TF) {
    TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
    TM[cbind(seq_len(length(TF)),TF)]=1 # factor indicator matrix  
    TM
    #NM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(NM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
  }
  # ******************************************************************************************************************************************************** 
  # Main body of rows design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************     
  hcf=HCF(replicates)
  tol=1.000001
  if (!is.data.frame(treatments)) {
    TF=data.frame(Treatments=as.factor(rep(rep(rep(1:sum(treatments)),rep(replicates,treatments)/hcf),hcf)))  
    model=" ~ Treatments"
    } else TF=treatments[rep(seq_len(nrow(treatments)), replicates), ,drop=FALSE]
  if (is.null(model)) model=paste("~",paste0(colnames(TF), collapse="*"))
  TF=data.frame(do.call(rbind,lapply(1:hcf,function(i) {TF[ sample( (1+(i-1)*nrow(TF)/hcf):(i*nrow(TF)/hcf) ), ,drop=FALSE]}))) 
  fnames=colnames(TF)
  rownames(TF) = seq(length=nrow(TF))
  nunits=nrow(TF)    
  if (length(columns)==0) columns=rep(1,length(rows))
  cumcols=c(1,cumprod(columns))
  cumrows=c(1,cumprod(rows))
  cumblocks=c(1,cumprod(rows*columns))
  strata=length(rows)
  if (is.null(searches)) searches=1+10000%/%nunits
  set.seed(seed)
  isrowcol=max(columns)>1
  blocksizes=nunits
  for (i in seq_len(strata)) 
    blocksizes=Sizes(blocksizes,i) 
  if (max(blocksizes) > min(blocksizes)) regBlocks=FALSE else regBlocks=TRUE
  blkdesign=data.frame(lapply(1:(strata+1),function(i) {gl( cumblocks[i], (cumblocks[strata+1]/cumblocks[i]), labels=unlist(lapply(1:cumblocks[i], function(j) {paste0("Block_",j)})))}))
  colnames(blkdesign)=unlist(lapply(0:(ncol(blkdesign)-1), function(j) {paste0("Stratum_",j,":Blocks")}))
  Stratum_0=as.factor(rep(1,nunits))
  Plots=as.factor(seq_len(nunits))
  if (!isrowcol) {
    rowdesign=data.frame(lapply(1:strata,function(i) {gl( (cumrows[i+1]), (cumblocks[strata+1]/cumblocks[i+1]), labels=unlist(lapply(1:cumrows[i+1], function(j) {paste0("Block_",j)})))}))
    colnames(rowdesign)=unlist(lapply(1:ncol(rowdesign), function(j) {paste0("Stratum_",j)}))
  } else if (isrowcol) {
    rowdesign=data.frame(lapply(1:strata,function(i) {gl( (rows[i]), (cumblocks[strata+1]/cumblocks[i]/rows[i]),cumblocks[strata+1],
                                                             labels=unlist(lapply(1:rows[i], function(j) {paste0("Row_",j)}))) }))
    coldesign=data.frame(lapply(1:strata,function(i) {gl( columns[i], (cumblocks[strata+1]/cumblocks[i+1]), cumblocks[strata+1],
                                                             labels=unlist(lapply(1:columns[i], function(j) {paste0("Col_",j)}))) }))
    rowdesign=data.frame(lapply(1:ncol(rowdesign), function(i){ interaction(blkdesign[,i], rowdesign[,i], sep = ":", lex.order = TRUE) }))
    coldesign=data.frame(lapply(1:ncol(coldesign), function(i){ interaction(blkdesign[,i], coldesign[,i], sep = ":", lex.order = TRUE) }))
    colnames(rowdesign)=unlist(lapply(1:ncol(rowdesign), function(j) {paste0("Stratum_",j,":Rows")}))
    colnames(coldesign)=unlist(lapply(1:ncol(coldesign), function(j) {paste0("Stratum_",j,":Cols")}))
    colDesign=data.frame(coldesign[rep(seq_len(length(blocksizes)),  blocksizes ),])  
    colDesign=data.frame(Stratum_0, colDesign)
  }
  blkDesign=data.frame(blkdesign[rep(seq_len(length(blocksizes)),  blocksizes ),])  
  rowDesign=data.frame(rowdesign[rep(seq_len(length(blocksizes)),  blocksizes ),])  
  rowDesign=data.frame(Stratum_0, rowDesign,Plots)
  regReps=identical(length(replicates),as.integer(1))
  orthoMain=(regReps && (replicates[1]==rows[1]))
    for ( i in 1:strata) {
      if (!isrowcol && rows[i]>1) TF=rowsOpt(TF,rowDesign[,i],rowDesign[,i+1])
      if ( isrowcol && rows[i]>1) TF=rowsOpt(TF,blkDesign[,i],rowDesign[,i+1],weighted,TM)
      if ( isrowcol && cols[i]>1) TF=colsOpt(TF,Model,blkDesign[,i],rowDesign[,i+1],colDesign[,i+1],blkDesign[,i+1],weighted,TM)
    }
  colnames(TF)=fnames
  if (  (max(replicates)>1 | is.data.frame(treatments))   && !isrowcol ) {
    rDesign=data.frame(do.call(cbind,lapply(1:ncol(rowDesign), function(r){ sample(nlevels(rowDesign[,r]))[rowDesign[,r]] }))) # Randomize
    rDesign=data.frame(rDesign,TF) 
    rDesign=rDesign[do.call(order, rDesign), ] # re-order
    blocksizes=table(rDesign[,ncol(rDesign)-ncol(TF)-1])[unique(rDesign[,ncol(rDesign)-ncol(TF)-1])]
    Design  = data.frame(rowdesign[rep(seq_len(length(blocksizes)),  blocksizes ),],seq_len(nunits), rDesign[,c((ncol(rDesign)-ncol(TF)+1):ncol(rDesign))])  # rebuild factor levels in order
    rownames(Design)=NULL
    if (!is.data.frame(treatments)) colnames(Design)=c(colnames(rowdesign),"Plots","Treatments") else
      colnames(Design)=c(colnames(rowdesign),"Plots",fnames)
    V=split(Design[,c((ncol(Design)-ncol(TF)+1):ncol(Design))],Design[,(ncol(Design)-ncol(TF)-1)])
    PlotsInBlocks=rep("",length(V))
    if (is.data.frame(treatments)) Plan=do.call(rbind, V)
    else Plan=as.data.frame(cbind(rowdesign,PlotsInBlocks, do.call(rbind, lapply(V, function(x){ length(x) =max(blocksizes); x }))))
    rownames(Plan)=NULL
  }
  if (max(replicates)>1 && isrowcol ) {
    rdf = data.frame(rowdesign,coldesign,blkdesign[,c(2:ncol(blkdesign))])[c(rbind(1:strata,(strata+1):(2*strata),(2*strata+1):(3*strata)))]
    rDesign=data.frame(rdf[rep(seq_len(length(blocksizes)), blocksizes),],Plots)
    rDesign=data.frame(lapply(1:ncol(rDesign), function(r){ sample(nlevels(rDesign[,r]))[rDesign[,r]]})) # Randomize
    rDesign=data.frame(rDesign, TF)
    rDesign=rDesign[do.call(order, rDesign), ] # re-order
    blocksizes=table(rDesign[,(ncol(rDesign)-2)])[unique(rDesign[,(ncol(rDesign)-2)])]
    Design  = data.frame(rdf[rep(seq_len(length(blocksizes)), blocksizes),],as.factor(rDesign[,ncol(rDesign)])) # rebuild factor levels in order
    
    colnames(Design)[ncol(Design)]="Treatments"
    rownames(Design)=NULL
    if (max(blocksizes)>1) {
      V=split(Design[,ncol(Design)],Design[,(ncol(Design)-1)])
      Plots_In_Blocks=rep("",length(V))
      Plan=as.data.frame(cbind(rdf[,-c(1:strata)*3],Plots_In_Blocks, do.call(rbind, lapply(V, function(x){ length(x) =max(blocksizes); x }))))
    } else {
      V=split(Design[,ncol(Design)],Design[,(ncol(Design)-3)] )
      plan=rdf[c(1:rows[strata])*columns[strata],-c( ncol(rdf)-1, ncol(rdf)),drop=FALSE]
      rownames(plan)=NULL
      Columns=rep("",nrow(plan))
      Plan=as.data.frame(cbind(plan,Columns, do.call(rbind, lapply(V, function(x){ length(x) =columns[strata]; x }))))
    }
    rownames(Plan)=NULL
  }

  # efficiencies
  if (!is.data.frame(treatments) && isrowcol) Efficiencies=RowColEfficiencies(Design)
  else if (!is.data.frame(treatments) && (max(replicates)==1)) Efficiencies=UnrepEfficiencies()
  else if (!is.data.frame(treatments)) Efficiencies=BlockEfficiencies(Design)
  else if (is.data.frame(treatments) && isrowcol) Efficiencies=BlockEfficiencies(Design)
  else if (is.data.frame(treatments) && !isrowcol) Efficiencies=factEffics(Design)
  row.names(Efficiencies)=NULL

  # omit single level row or column strata in row and column designs
  if (isrowcol && !is.data.frame(treatments) ) {
    Design=Design[,-c((1:strata)*3)]
    Design[c(which(as.numeric(rbind(rows,columns))==1))]= list(NULL) 
    Plan[c(which(as.numeric(rbind(rows,columns))==1 ))]= list(NULL) 
  } 
  # treatment replications
  if (!is.data.frame(treatments) ) {
    TreatmentsTable=data.frame(table(Design[,ncol(Design)]))
    TreatmentsTable[]=lapply(TreatmentsTable, as.factor) 
    colnames(TreatmentsTable)=c("Treatments","Replicates")
  } else {
    TreatmentsTable=NULL
    Plan=NULL
  }
  list(Treatments=TreatmentsTable,Efficiencies=Efficiencies,Plan=Plan,Design=Design,Seed=seed,Searches=searches,Jumps=jumps) 
} 