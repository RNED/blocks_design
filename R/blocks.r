#' @title Block designs 
#' 
#' @description
#' 
#' Constructs nested block designs for unstructured treatments with arbitrary replication, 
#' not necessarily all equal, and arbitrary depth of nesting.
#' 
#' @details
#' 
#' The algorithm optimizes the allocation of treatments to blocks for nested blocks designs.
#' Treatments can have any arbitrary replication not necessarily all equal and the blocks can be simple main blocks or can
#' be nested blocks with any feasible depth of nesting. Block sizes in any particular stratum are always equal or differ, at most, by a single unit. 
#' 
#' Treatments are defined by a partition into treatment sets where all treatments in the
#' same set have the same replication. The required partition is defined in the \code{treatments} list and the 
#' corresponding replication for each treatment set is defined in the \code{replicates} list. The two lists must be the same length 
#' and must be in matching order. Treatments are numbered consecutively according to the ordering of the treatment sets but 
#' different sets with the same replication can be used if arbitrary numbering is required. Single replicate treatments sets are permitted provided 
#' that not every treatment in the design is unreplicated.
#' 
#' The blocks design is defined by the \code{blocklevels} list, which contains one block level for each stratum of the design. 
#' The first level is the number of main blocks and the successive levels, if any, are the numbers of nested sub-blocks in each succesive strata of
#' the design. The length of the list is the total number of strata in the design. \code{blocklevels} has a single default value
#' which is the highest common factor of the replication levels and which gives a main blocks design with a maximal set of orthogonal main blocks.  
#' 
#' The blocks design is optimized hierarchically with the top level blocks optimized unconditionally and the succesive nested blocks, if any, 
#' optimized hierarchically from top-down with the blocks of each new nested stratum optimized within the blocks of each preceding stratum.
#' Thus nested blocks are optimized conditionally with respect to any higher level blocks in the design but unconditionally with respect to any lower level blocks.
#'  
#' If the blocks in the top stratum have k replicates with v**2 equally replicated treatments in blocks of size v
#' and k <= 3 for any v, or k <= v+1 for prime or prime-power v, or k <= 4 for v = 10, 
#' they are regular lattice blocks and are constructed algebraically. All other blocks are constructed by a D-optimal swapping algorithm. 
#'  
#' The algorithm searches for a local maxima in each stratum by making improving swaps between blocks nested within existing blocks 
#' until no further improving swaps can be found. If the \code{searches} parameter is greater than one the algorithm 
#' then escapes the current local maxima by making one or more random jumps according to the \code{jumps} parameter and then finds another local maxima. 
#' The process continues for the required number of searches and then the best overall design is returned. 
#' 
#' The design outputs include the Design showing the allocation of blocks and treatments to plots, a Plan showing a schematic 
#' layout of treatments in blocks, a set of Incidence matrices showing the incidences of treatments and blocks in each stratum, 
#' an Efficiencies table showing the final achieved A- and D-efficiencies of the overall best design together with A-upper bounds, where available and 
#' a Searches_Log showing the progress of the search algorithm for design improvements in each stratum of the design. 
#' The Searches_Log shows the number of searches, the D-efficiency and the A-efficiency for each improving swap in the optimization process
#'  and can be useful for assessing whether further searches are needed or whether the \code{jumps} parameter setting
#'  is appropritae for a particular design search.
#'  
#'  Designs are fully randomised with small blocks randomised within large blocks and treatment plots randomised within blocks.  
#'  
#' The \code{searches} parameter is the number of local optima searched during an optimization. 
#' Increasing the number of searches may improve the efficiency of a design but
#'  will also increase the search time.
#'  
#' The \code{seed} parameter is an integer used to initialize the random number generator. The 
#'  default is a random integer but any fixed positive integer can be used instead,if required.   
#' 
#' Blocks and treatments are fully randomized within the constraints of a nested blocks design.
#' 
#' @param treatments a list giving a partition of the total number of treatments into 
#' sets where all treatments in the same set have the same replication.   
#' 
#' @param replicates a list assigning the replication level of each set in the \code{treatments} list. 
#' 
#' @param blocklevels a list of block levels where the first level is the number of main blocks and the remaining
#' levels, if any, are the succesive nested levels of a hierarchy of nested sub-blocks.
#' The default is an orthogonal main blocks design.
#' 
#' @param seed an optional integer parameter for initializing the random number generator. The default 
#' is a random integer seed in the range 1:100000.
#' 
#' @param searches an optional integer parameter for the number of local optima searched during an optimization. 
#' The default is the minimum of 32 or the ceiling of 4096 divided by the number of units.
#' 
#' @param jumps an optional integer parameter for the number of random jumps used to escape a local maxima in each stratum. 
#' The default is a single jump.
#' 
#' @return  
#' \item{Design}{Data frame showing the listing of treatments allocated to blocks}
#' \item{Plan}{Data frame showing a plan of treatments allocated to sub-plots within blocks}
#' \item{Incidences}{Blocks-by-treatments incidence matrices, one for each stratum of the design}
#' \item{Efficiencies}{Data frame showing the achieved efficiencies for each stratum of the design together with an A-efficiency upper-bound, where available}
#' \item{Searches_Log}{Data frames showing the number of searches for each progressive improvement in design efficiency for each stratum of the design} 
#' \item{seed}{Numerical seed for random number generator}
#' \item{searches}{Maximum number of searches in each stratum}
#' \item{jumps}{Number of jumps to escape a local maxima in each stratum}
#'
#' @references
#' 
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. http://CRAN.R-project.org/package=crossdes
#' 
#' @examples
#' 
#' # 3 treatments with 2 reps, 2 treatments with 4 reps, 4 treatments with 3 reps 
#' # the hcf of the replication numbers is 1 and the default design is a completely randomized design 
#' blocks(treatments=c(3,2,4),replicates=c(2,4,3))
#' 
#' # 50 treatments with 4 reps in 4 complete randomized blocks 
#' blocks(treatments=50,replicates=4)
#' 
#' # as above but with 4 main blocks and 5 nested blocks within each main block 
#' blocks(treatments=50,replicates=4,blocklevels=c(4,5))
#' 
#' # as above but with 20 additional single replicate treatments, one to each block
#' blocks(treatments=c(50,20),replicates=c(4,1),blocklevels=c(4,5))
#' 
#' # 64 treatments with 2 reps and 2 main blocks with five 2-level nested factors   
#'  blocks(treatments=64,replicates=2,blocklevels=c(2,2,2,2,2,2),searches=4)
#' 
#' # concurrence matrices of 36 treatments with 3 reps and 3 main blocks with 6 nested blocks
#' crossprod(blocks(treatments=36,replicates=3,blocklevels=c(3,6))$Incidences[[2]])
#' 
#' # concurrence matrix for 13 treatments with 4 reps and 13 treatments with one rep in 13 blocks 
#' crossprod(blocks(c(13,13),c(4,1),13,searches=100)$Incidences[[1]])
#' 
#'          
#' @export
#' 
blocks = function(treatments, replicates, blocklevels=NULL, searches=NULL, seed=NULL,jumps=NULL) { 

  #******************************************************** sizes ************************************************************************************ 
  Sizes=function(sizes,blocklevels) { 
    for  (i in 1:length(blocklevels)) {    
      newsizes=NULL
      for (z in 1: length(sizes)) 
        newsizes=c(newsizes, rep(sizes[z] %/% blocklevels[i], blocklevels[i]) + c( rep(1, sizes[z] %% blocklevels[i]), rep(0,(blocklevels[i]-sizes[z] %% blocklevels[i])))) 
      sizes=newsizes
    }   
    sizes 
  } 
     
  #******************************************************** Primality test ************************************************************************************
  isPrime=function(v) {
    if (v <= 3)  return(TRUE)
    else if (v %% 2 == 0 | v %% 3 == 0) return(FALSE)
    else if (v<25) return(TRUE)
    else 
      for(i in  6*rep(1:floor((sqrt(v)+1)/6)) )
        if( v %% (i-1) == 0 | v %% (i+1) == 0) return(FALSE) 
    return(TRUE)
  }      
  
  #******************************************************** Block contrasts ************************************************************************************
  BlockContrasts=function(MF,BF) {
    BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
    BM[cbind(1:length(BF),as.numeric(BF))]=1 # factor indicator matrix
    BM=BM[,rep(c(rep(TRUE,((nlevels(BF)/nlevels(MF))-1)),FALSE),nlevels(MF)),drop=FALSE]
    for (i in 1: nlevels(MF)) 
      BM[MF==i,]=scale(BM[MF==i,], center = TRUE, scale = FALSE) # scaling within main blocks
    BM
  } 
  
  #******************************************************** Treatment contrasts**********************************************************************************
  TreatContrasts=function(MF,TF) {
    TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
    TM[cbind(1:length(TF),as.numeric(TF))]=1 # factor indicator matrix	
    TM=TM[,c(rep(TRUE,(nlevels(TF)-1)),FALSE),drop=FALSE]	
    for (i in 1:nlevels(MF)) 
      TM[MF==i,]=scale(TM[MF==i,] , center = TRUE, scale = FALSE)
    TM
  }
  
  #******************************************************** Updates variance matrix ************************************************************************************
  UpDate=function(M11,M22,M12,ti,tj,bi,bj,TF,BF) {  
    m11=M11[ti,ti]+M11[tj,tj]-M11[tj,ti]-M11[ti,tj]
    m22=M22[bi,bi]+M22[bj,bj]-M22[bi,bj]-M22[bj,bi]
    m12=M12[ti,bi]-M12[tj,bi]-M12[ti,bj]+M12[tj,bj]   
    f = sqrt(2+m11+m22-2*m12)
    m = f/sqrt(1-2*m12-m11*m22+m12*m12)/2 
    Z1 = (M12[,bi]-M12[,bj]-M11[,ti]+M11[,tj])/f     
    Z2 = (M22[bi,]-M22[bj,]-M12[ti,]+M12[tj,])/f 
    W1 = (M11[,ti]-M11[,tj]+M12[,bi]-M12[,bj] - Z1*(m22-m11)/f)*m
    W2 = (M12[ti,]-M12[tj,]+M22[bi,]-M22[bj,] - Z2*(m22-m11)/f)*m
    M11 = M11 - tcrossprod(Z1) + tcrossprod(W1)
    M22 = M22 - tcrossprod(Z2) + tcrossprod(W2)
    M12 = M12 - tcrossprod(Z1,Z2) + tcrossprod(W1,W2) 
    list(M11=M11,M22=M22,M12=M12)
  }   
  
  #********************************************************Determinants of jumps using samples of increasing size***********************************************
  D_Max=function(M11,M22,M12,TF,MF,BF) {      
    locrelD=1
    mainSizes=tabulate(MF)
    nSamp=pmin(rep(8,nlevels(MF)),mainSizes)
    mainBlocks=split(rep(1:length(TF)),MF)  
    repeat {
      improved=FALSE
      for (k in 1:nlevels(MF)) {
        S=sort(sample(mainBlocks[[k]],nSamp[k]))  
        TT=2*M11[TF[S],TF[S],drop=FALSE]-tcrossprod(M11[cbind(TF[S],TF[S])]+rep(1,nSamp[k])) + tcrossprod(M11[cbind(TF[S],TF[S])]) + 1
        BB=2*M22[BF[S],BF[S],drop=FALSE]-tcrossprod(M22[cbind(BF[S],BF[S])]+rep(1,nSamp[k])) + tcrossprod(M22[cbind(BF[S],BF[S])]) + 1
        TB=M12[TF[S],BF[S],drop=FALSE]+t(M12[TF[S],BF[S],drop=FALSE])-tcrossprod(M12[cbind(TF[S],BF[S])]+rep(1,nSamp[k]))+tcrossprod(M12[cbind(TF[S],BF[S])]) + 2
        dMat=TB**2-TT*BB
        sampn=which.max(dMat)   
        i=1+(sampn-1)%%nSamp[k]
        j=1+(sampn-1)%/%nSamp[k]
        relD=dMat[i,j]
        if ( !isTRUE(all.equal(relD,1)) & relD>1) {
          improved=TRUE
          locrelD=locrelD*relD
          up=UpDate(M11,M22,M12,TF[S[i]],TF[S[j]], BF[S[i]], BF[S[j]], TF,BF)
          M11=up$M11
          M22=up$M22
          M12=up$M12
          TF[c(S[i],S[j])]=TF[c(S[j],S[i])]
        }
      } 
      if (improved) next
      if (sum(nSamp) < min(length(TF),512))
          nSamp=pmin(mainSizes,2*nSamp)
       else 
         break
    }  
    list(M11=M11,M22=M22,M12=M12,TF=TF,locrelD=locrelD)
  }
  
  #**************************** Calculates A-optimality *******************************************************
  optEffics=function(TF,BF,ntrts,nblks)   {
    NN= t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF)))  
    if (ntrts<=nblks) 
      e=eigen( (diag(ntrts)-crossprod(NN)), symmetric=TRUE, only.values = TRUE)$values[1:(ntrts-1)]     
    else       
      e=c(rep(1,(ntrts-nblks)),eigen((diag(nblks)-tcrossprod(NN)), symmetric=TRUE, only.values = TRUE)$values[1:(nblks-1)])    
    aeff = 1/mean(1/e) 
    deff = exp(sum(log(e))/(ntrts-1))
    c(deff,aeff)
  }
  
  #**************************** General optimization using annealing with multiple searches *******************************************************
  Optimise=function(TF,BF,MF,M11,M22,M12,searches,Iterations,jumps)  {
    relD=1
    globrelD=0
    globTF=TF
    treps=as.integer(tabulate(TF))
    breps=as.integer(tabulate(BF))
    bound=NA
    if (identical(max(treps),min(treps)) & identical(max(breps),min(breps))  )
        bound=upper_bounds(length(TF),nlevels(TF),nlevels(BF)) 
    for (r in 1 : searches) {
      dmax=D_Max(M11,M22,M12,TF,MF,BF)  
      relD=relD*dmax$locrelD
      TF=dmax$TF
      M11=dmax$M11
      M22=dmax$M22
      M12=dmax$M12  
        if (relD>globrelD) {
        globTF=TF
        globrelD=relD
        effics=optEffics(globTF,BF,nlevels(TF),nlevels(BF)) 
        Iterations=rbind(Iterations,c(r,effics))
        if (isTRUE( all.equal(bound,effics[2]))) break
        }
      if (r==searches) break
      for (iswap in 1 : jumps) {   
        for (icount in 1 : 100) {
          s=sample(rep(1:length(TF))[MF==sample(nlevels(MF),1)],2)
          if ( identical(TF[s[1]],TF[s[2]]) | identical(BF[s[1]],BF[s[2]])  ) next
          dswap = (1+M12[TF[s[1]],BF[s[2]]]+M12[TF[s[2]],BF[s[1]]]-M12[TF[s[1]],BF[s[1]]]-M12[TF[s[2]],BF[s[2]]])**2-
            (2*M11[TF[s[1]],TF[s[2]]]-M11[TF[s[1]],TF[s[1]]]-M11[TF[s[2]],TF[s[2]]])*(2*M22[BF[s[1]],BF[s[2]]]-M22[BF[s[1]],BF[s[1]]]-M22[BF[s[2]],BF[s[2]]])  
          if (!isTRUE(all.equal(dswap,0)) & dswap>0) {
            relD=relD*dswap
            up=UpDate(M11,M22,M12,TF[s[1]],TF[s[2]], BF[s[1]], BF[s[2]],TF,BF)
            M11=up$M11
            M22=up$M22
            M12=up$M12
            TF[c(s[1],s[2])]=TF[c(s[2],s[1])]  
            break
          }
        }
      }
    } 
    list(TF=globTF,Iterations=Iterations)
  } 
  
  #******************************************************** Initializes design***************************************************************************************
  GenOpt=function(TF,BF,MF,searches,Iterations,jumps) {   
    TB=TreatContrasts(MF,TF)
    NB=BlockContrasts(MF,BF) 
    DD=crossprod(cbind(TB,NB))
    count=0
      while ( !identical(qr(DD)$rank,ncol(DD)) & count<100 ) 
      {
        rand=sample(1:length(TF))
        TF=TF[rand][order(MF[rand])] 
        TB=TB[rand,][order(MF[rand]),]
        DD=crossprod(cbind(TB , NB))
        count=count+1
      }
    if (count==100) return(TF)   
    V=chol2inv(chol(DD))
    M11=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))	
    M22=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
    M12=matrix(0,nrow=nlevels(TF),ncol=nlevels(BF))
    M11[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[1:(nlevels(TF)-1),1:(nlevels(TF)-1),drop=FALSE]
    M12[1:(nlevels(TF)-1),1:(ncol(V)-nlevels(TF)+1)]=V[1:(nlevels(TF)-1),nlevels(TF):ncol(V),drop=FALSE]
    M22[1:(ncol(V)-nlevels(TF)+1),1:(ncol(V)-nlevels(TF)+1)]=V[nlevels(TF):ncol(V),nlevels(TF):ncol(V),drop=FALSE]
    perm=order(order((1:nlevels(BF))%%(nlevels(BF)/nlevels(MF))==0)   )
    M12=M12[,perm]
    M22=M22[perm,perm]	
    opt=Optimise(TF,BF,MF,M11,M22,M12,searches,Iterations,jumps)
    list(TF=opt$TF,Iterations=opt$Iterations)
  }
  
  #******************************************************** HCF of replicates************************************************************************************
  HCF=function(replevs)  {
    replevs=sort(replevs)
    v=(c(replevs[1],NULL))
    if (length(replevs)>1) 
      for (i in 2: length(replevs)) {
        v[2]=replevs[i] 
        while (v[2]%%v[1] != 0) v = c(v[2]%%v[1], v[1]) 
      }
    v[1]
  }
 
  #********************************************************  algorithm   ************************************************************************************
    optTF=function(Design,treatlevs,replevs,searches,jumps)  {
    nunits=nrow(Design)
    ntrts=sum(treatlevs)
    hcf=HCF(replevs)
    v=sqrt(ntrts)
    strata=ncol(Design)-2
    Iterations=vector("list", strata)
    for (i in 1:strata)
      Iterations[[i]]=matrix(0,nrow=1,ncol=3) 
    orthbsize=nunits/hcf 
    ortho=0
    for (i in 1 : strata) 
      if (all( tabulate(Design[,i+1]) %% orthbsize == 0)) ortho=i else  break
    reglat=( all(replevs==replevs[1]) & max(Design[,i+1])==v*replevs[1] & identical(v,ntrts%/%v) )
    pp_trts=c(16,64,256,1024,4096,16384,81,729,6561,625,2401)  
    simplelattice = (reglat &  replevs[1]<4 )
    primelattice =  (reglat &  replevs[1]<(v+2)  & isPrime(v))
    ppowerlattice= (reglat  &  replevs[1]<(v+2)  &  ntrts%in% pp_trts)
    lattice100 =(reglat & v==10  & replevs[1]<5 )  
    # treps is the vector of treatment replications for the minimum orthogonal block size
    treps=rep(replevs,treatlevs)/hcf  
    TF=rep(rep(1:ntrts,treps),hcf)
    if (ortho<strata) {
      for (i in (ortho+1) : strata) { 
        if ( i==(ortho+1)  & simplelattice) {		
          TF=c(rep(1:(v*v)), rep(1:(v*v))[order(rep(0:(v-1),v))])
          if (replevs[1]>2) {
            set=NULL
            for (j in 0: (v-1)) 
              for (k in 0: (v-1)) 
                set=c(set, (j+k)%%v )
            TF=c(TF, rep(1:(v*v))[order(set)])
          }		
        } else if ( i==(ortho+1)  & primelattice ) {  		
          TF=c(rep(1:(v*v)), rep(1:(v*v))[order(rep(0:(v-1),v))])
          for (z in 1: (replevs[1]-2)) {
            set=NULL
            for (j in 0: (v-1)) 
              for (k in 0: (v-1)) 
                set=c(set,(j+k*z)%%v)
            TF=c(TF, rep(1:(v*v))[order(set)])		
          }	
        } else if ( i==(ortho+1) & ppowerlattice ) {	
          prime= c(2,2,2,2,2,2,   3,3,3,  5,7)[which(pp_trts==ntrts)]
          ppower=c(2,3,4,5,6,7,   2,3,4,  2,2)[which(pp_trts==ntrts)]
          mols=crossdes::MOLS(prime,ppower)			
          TF=c(rep(1:(v*v)), rep(1:(v*v))[order(rep(0:(v-1),v))])
          for (i in 1: (replevs[1]-2))
            TF=c(TF, rep(1:(v*v))[order(    as.numeric(mols[,,i]) ) ])
        } else if (i==(ortho+1) &lattice100) {
          TF=c(rep(1:(v*v)), rep(1:(v*v))[order(rep(0:(v-1),v))])
          if (replevs[1]>2)  { 
            tens1=c(
              1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
              5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8)
            TF=c(TF, rep(1:(v*v))[order(tens1)]) 
          }
          if (replevs[1]==4) {
            tens2=c(
              1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
              6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6)
            TF=c(TF, rep(1:(v*v))[order(tens2)])
          }
        } else {	
          opt=GenOpt(as.factor(TF),as.factor(Design[,(i+1)]),as.factor(Design[,i]),searches,Iterations[[i]],jumps)
          TF=opt$TF
          Iterations[[i]]=opt$Iterations
        }
      }
    }
    list(TF=TF,Iterations=Iterations)  
  }
  
  #******************************************************** replaces single rep treatments *************************************************************************** 
  fullDesign=function(Design,treatments,replicates,blocksizes,blocklevels) {
    strata=ncol(Design)-3
    TF=as.factor(Design[,ncol(Design)])
    nunits=sum(treatments*replicates)
    ntrts=sum(treatments)
    redtrts=nlevels(TF)
    TF=as.factor(c(TF,((redtrts+1):ntrts)[sample(ntrts-redtrts)])) 
    trtlabs=NULL  
    extlabs=NULL
    index=0
    for (i in 1 : length(treatments)) {
      if (replicates[i]>1) 
        trtlabs=c(trtlabs,  (index+1):(index+treatments[i]) )
      else extlabs=c(extlabs,  (index+1):(index+treatments[i]) )
        index=index+treatments[i] 
    }
    trtlabs=c(trtlabs,extlabs)
    levels(TF)=trtlabs
    TF=as.numeric(levels(TF))[TF]
    BF=c( rep( 1:length(blocksizes),blocksizes))
    newblocksizes=Sizes(nunits,blocklevels)
    BF=c(BF,  rep( 1:length(blocksizes),(newblocksizes-blocksizes) ) )
    # full TF in blocks
    TF=TF[order(BF)]
    Design = matrix(1,nrow=nunits,ncol=(strata+3))
    for (r in 1 : strata) 
      Design[,r+1]=rep(facMat[,r],newblocksizes)
    Design[,r+2]=rep(1:nunits)
    Design[,r+3]=as.factor(TF)  
    Design
  } 
  
  #******************************************************** A-efficiencies ************************************************************************************
  A_Efficiencies=function(Design)  {
    strata=ncol(Design)-2
    treps=tabulate(Design$Treatments)
    effics=matrix(1,nrow=strata,ncol=2)
    aeff=rep(1,strata) 
    deff=rep(1,strata)  
    bounds=rep(NA,strata) 
    blocks=rep(0,strata)  
    for (i in 1:strata) { 
      blocks[i]=nlevels(Design[,i])
      breps=tabulate(Design[,i])
      if ( all(treps==treps[1]) & all(breps==breps[1]))
        bounds[i]=upper_bounds(nrow(Design),nlevels(Design$Treatments),blocks[i])    
      if (nlevels(Design$Treatments)>1 & nlevels(Design[,i])>1)
        effics[i,]=optEffics(Design$Treatments,Design[,i],nlevels(Design$Treatments),blocks[i])  
    }
    efficiencies=as.data.frame(cbind(blocks, effics, bounds))    
    colnames(efficiencies)=c("Blocks","D-Efficiencies","A-Efficiencies", "A-Upper Bounds")
    rnames=c("Main")
    if (strata>1)
      for (i in 1 : (strata-1)) rnames=c(rnames,paste("Sub",i))
    rownames(efficiencies)=rnames 
    efficiencies
  }
   
  #******************************************************** Randomizes blocks within strata************************************************************************************ 
  randBlocks=function(Design,facMat) {
    for (r in 1 : (ncol(Design)-1) )
      Design[,r]=as.numeric( sample(nlevels( Design[,r]) ))[Design[,r]]  
    Design=Design[ do.call(order, Design), ] 
    blocksizes=tabulate(as.numeric( order(unique(Design[,ncol(Design)-2])))[Design[,ncol(Design)-2]])
    for (r in 1 : (ncol(Design)-2) ) 
      Design[,r]=rep(facMat[,r],blocksizes)
    Design[,(ncol(Design)-1)]=rep(1:nrow(Design))
    Design[]=lapply(Design, factor) 
    Design
  }
 
  #******************************************************** Plan output************************************************************************************
  Plan=function(Design,facMat,stratumnames)  {
    strata=ncol(Design)-2
    bSizes=c(0,tabulate(Design[,strata]))
    nblocks=length(bSizes)-1
    plotTrts=matrix(nrow=nblocks,ncol=max(bSizes)) 
    for (i in 1:nblocks) 
      plotTrts[i, (1 : bSizes[i+1])] = Design[(1 + sum(bSizes[1:i])) : sum(bSizes[1:(i+1)]) , strata+2]  
    plotTrts[is.na(plotTrts)]  = " "
    Plan=as.data.frame(cbind(facMat, rep(" ",nblocks), plotTrts))
    stratumnames=c(stratumnames, "Sub_plots", 1:ncol(plotTrts) )
    colnames(Plan)=stratumnames
    Plan
  }
  
 #******************************************************** Validates inputs************************************************************************************
 testInputs=function(treatments,replicates,blocklevels,searches,seed,jumps) {  
   if (missing(treatments) | missing(replicates) )  
     return(" Treatments or replicates not defined ")   
   if (is.null(treatments) | is.null(replicates))  
     return(" Treatments or replicates list is empty ")   
   if (anyNA(treatments) | anyNA(replicates) ) 
     return(" NA values not allowed")
   if (!all(is.finite(treatments)) | !all(is.finite(replicates)) | !all(!is.nan(treatments)) | !all(!is.nan(replicates))) 
     return(" Treatments and replicates can contain only finite integers ")
   if ( length(treatments)!=length(replicates) ) 
     return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))
  if (!all(treatments>=1)) 
     return("Treatments must be integers greater than zero")
  if (!all(replicates>=1)) 
    return("Replicates must be integers greater than zero")  
  if (all(replicates==1)) 
    return("Not all treatments can have only a single replication" )  
   if (!is.null(blocklevels)) {
     if (anyNA(blocklevels) ) return(" NA blocklevels values not allowed") 
     if (!all(is.finite(blocklevels)) | !all(!is.nan(blocklevels)) ) return(" Blocklevels can contain only finite integers ")
     if (min(blocklevels)<1) return (" Blocklevels must be at least one ")
   }
   if (!is.null(searches)) {
     if (anyNA(searches) ) return(" NA searches values not allowed") 
     if ( !all(is.finite(searches)) | !all(!is.nan(searches))) return(" Searches must be a finite integer ") 
     if (searches<1)  return(" Repeats must be at least one ")   
   }  
  
  if (!is.null(jumps)) {
    if (anyNA(jumps) ) return(" NA jumps values not allowed") 
    if ( !all(is.finite(jumps)) | !all(!is.nan(jumps))) return(" jumps must be a finite integer ") 
    if (jumps<1)  return(" Random jumps must be at least one ")   
  }  
  
   if (!is.null(seed)) {
     if (anyNA(seed) ) return(" NA seed values not allowed") 
     if ( !all(is.finite(searches)) | !all(!is.nan(searches))) return(" Seed must be a finite integer ") 
     if (seed<1)  return(" Seed must be at least one ")   
   }  
   if (  sum(treatments*replicates) < (prod(blocklevels) + sum(treatments)) ) 
     return("Design cannot be fitted :  too many blocks and treatments for the available plots")  
   return(TRUE)
 }
 
 #********************************************************Iterationss design ************************************************************************************ 
 testout=testInputs(treatments,replicates,blocklevels,searches,seed,jumps) 
  if (!isTRUE(testout)) stop(testout)
  if (is.null(seed)) seed=sample(1:100000,1)
  set.seed(seed) 
 if (is.null(jumps)) jumps=1
  # omit any single replicate treatments here 
   treatlevs=treatments[replicates>1]
   replevs = replicates[replicates>1]
  if (is.null(blocklevels)) 
    blocklevels=HCF(replevs)
  nunits=sum(treatlevs*replevs) 
  if (is.null(searches)) 
   searches=min(64, ceiling(4096/nunits))
 if (!all(blocklevels==1))
    blocklevels=blocklevels[blocklevels>1]
 else
   blocklevels=1
  strata=length(blocklevels)
 
 for (i in 1:strata)
  blocksizes=Sizes(nunits,blocklevels)
 totblocks=prod(blocklevels)
 facMat= matrix(nrow=totblocks,ncol=strata)
  for (r in 1 : strata) 
    facMat[,r]=gl(prod(blocklevels[1:r]),totblocks/prod(blocklevels[1:r])  )  
 Design= matrix(1,nrow=nunits,ncol=(strata+2))
  for (r in 1 : strata) 
    Design[,r+1]=rep(facMat[,r],blocksizes)
 Design[,r+2]=rep(1:nunits)
 opt=optTF(Design,treatlevs,replevs,searches,jumps) 
 TF=opt$TF
 Iterations=opt$Iterations
  Design=cbind(Design,as.factor(TF)) 
  # add back single replicate treatments here 
  if (!all(replicates>1) )
   Design= fullDesign(Design,treatments,replicates,blocksizes,blocklevels) 
 Design=as.data.frame(Design)[,c(2:ncol(Design))] 
 Design[]=lapply(Design, factor)   
 # randomization
  Design=randBlocks(Design,facMat)
  stratumnames=c("Main_blocks")
  if (strata>1)
    for (i in 1:(strata-1))
      stratumnames=c(stratumnames,paste("Sub",i,"_blocks", sep=""))  
  colnames(Design)=c(stratumnames,"Sub-plots","Treatments")   
  rownames(Design) = NULL 
  # Incidence matrix for each stratum
  Incidences=vector(mode = "list", length =strata )
  for (i in 1:strata)
    Incidences[[i]]=table( Design[,i] ,Design[,strata+2])  
 names(Incidences)=stratumnames
 
 names(Iterations)=stratumnames
  plan=Plan(Design,facMat,stratumnames)
  efficiencies=A_Efficiencies(Design)
  
 for (i in 1:strata) {
   if (nrow(Iterations[[i]])==1)
     Iterations[[i]]=rbind(  Iterations[[i]] ,c(1,efficiencies[i,2],efficiencies[i,3] )  )
  Iterations[[i]] = Iterations[[i]] [  2:nrow(Iterations[[i]]),    ,drop=FALSE]
  Iterations[[i]]=as.data.frame(Iterations[[i]] )
  colnames(Iterations[[i]])=c("Searches","D-efficiency","A-efficiency")  
 }
 
  list(Design=Design,Plan=plan,Incidences=Incidences,Efficiencies=efficiencies,Progress_Log=Iterations,Seed=seed,Searches=searches,Jumps=jumps) 
} 