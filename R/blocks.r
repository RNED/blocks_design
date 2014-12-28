#' @title Block designs 
#' 
#' @description
#' 
#' Constructs nested block designs for unstructured treatments with arbitrary replication, 
#' not necessarily all equal, and arbitrary depth of nesting.
#' 
#' @details
#' 
#' \code{blocks(...)} provides a top-down optimization of a nested blocks design where the top stratum blocks are optimized unconditionally and
#'  the blocks of any nested strata are optimized conditionally within the blocks of each preceding stratum. 
#'  
#' If the blocks in the top stratum have k replicates with v**2 equally replicated treatments in blocks of size v
#' and k <= 3 for any v, or k <= v+1 for prime or prime-power v, or k <= 4 for v = 10, 
#' then the stratum is a lattice block design and is constructed algebraically.  
#'
#' All other blocks are constructed algoritmically by a D-optimality algorithm that makes 
#' improving swaps between nested blocks within containing blocks until no further improvement is possible.
#'  
#' \code{treatments} is a list of sets where the sum of the sets is the required number of treatments 
#' and the treatments in any one set are all equally replicated. 
#' 
#' \code{replicates} is a list of replication numbers for sets in the \code{treatments} list. 
#' Treatments are numbered consecutively according to the order of the sets
#' and treatments with the same replication can be split between two or more sets if required. 
#'  
#' \code{blocklevels} is a list of nested block levels for the succesive nested blocks strata of the design. 
#' The first level is the number of main blocks 
#' and the successive levels, if any, are the numbers of sub-blocks in the succesive strata of
#' the nested blocks design.
#' The length of the list is the number of strata  and the 
#' running products of the levels are the total blocks in each successive stratum of the
#' design. Blocks in the same stratum are always equal in size or differ by, at most, a
#' single unit. The default is the highest common factor of the replication levels, 
#' which gives a main blocks design with a maximal set of complete orthogonal main blocks. 
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
#' @param searches an optional integer for the number of local optima searched during an optimization. 
#' The default is the minimum of 32 or the ceiling of 4096 divided by the number of units.
#' 
#' @param seed an optional integer seed for initializing the random number generator. The default 
#' is a random seed.
#' 
#' @return  
#' \item{Design}{Data frame showing the listing of treatments allocated to blocks}
#' \item{Plan}{Data frame showing a plan of treatments allocated to sub-plots within blocks}
#' \item{Incidences}{List of blocks-by-treatments incidence matrices, one for each stratum of the design}
#' \item{Efficiencies}{Data frame showing the A-efficiency factor for each stratum of the design together with an upper bound, where available}
#' \item{seed}{Numerical seed for random number generator}
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
#' blocks(treatments=64,replicates=2,blocklevels=c(2,2,2,2,2,2),searches=4)
#' 
#' # concurrence matrices of 36 treatments with 3 reps and 3 main blocks with 6 nested blocks
#' crossprod(blocks(treatments=36,replicates=3,blocklevels=c(3,6))$Incidences[[2]])
#' 
#' # concurrence matrix for 13 treatments with 4 reps and 13 treatments with one rep in 13 blocks 
#' crossprod(blocks(c(13,13),c(4,1),13,searches=100)$Incidences[[1]])
#' 
#'     
#'       
#' @export
#' 

blocks = function(treatments, replicates, blocklevels=NULL, searches=NULL, seed=NULL) { 

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
  
  #********************************************************Determinants of swaps using samples of increasing size***********************************************
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
        if ( relD>1.000001 ) {
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
  
     
  #**************************** General optimization using annealing with multiple searches *******************************************************
  Optimise=function(TF,BF,MF,M11,M22,M12,searches,Iterations)  {
    relD=1
    globrelD=0
    treps=tabulate(TF)
    breps=tabulate(BF)
    if ( all(treps==treps[1]) )
      if (length(TF)%%nlevels(BF) == 0) 
        bound=upper_bounds(length(TF),nlevels(TF),nlevels(BF)) 
    else
      bound=NA
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
        test=optEffics(globTF,BF,treps,breps,nlevels(TF),nlevels(BF)) 
        Iterations=rbind(Iterations,c(r,test$deff,test$aeff))
        if (isTRUE( all.equal(bound,test$aeff))) break
        }
      if (r==searches) break
      for (iswap in 1 : 5) {
      for (icount in 1 : 100) {
          s=sample(rep(1:length(TF))[MF==sample(nlevels(MF),1)],2)
          if ( TF[s[1]]==TF[s[2]]| BF[s[1]]==BF[s[2]]) 	next
          dswap = (1+M12[TF[s[1]],BF[s[2]]]+M12[TF[s[2]],BF[s[1]]]-M12[TF[s[1]],BF[s[1]]]-M12[TF[s[2]],BF[s[2]]])**2-
            (2*M11[TF[s[1]],TF[s[2]]]-M11[TF[s[1]],TF[s[1]]]-M11[TF[s[2]],TF[s[2]]])*(2*M22[BF[s[1]],BF[s[2]]]-M22[BF[s[1]],BF[s[1]]]-M22[BF[s[2]],BF[s[2]]])  
          if ( isTRUE(all.equal(dswap,0)) | isTRUE(all.equal(dswap,1)) ) next
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
    list(TF=globTF,Iterations=Iterations)
  } 
  
  #******************************************************** Initializes design***************************************************************************************
  GenOpt=function(TF,BF,MF,searches,Iterations) {   
    TB=TreatContrasts(MF,TF)
    NB=BlockContrasts(MF,BF) 
    DD=crossprod(cbind(TB,NB))
    count=0
      while ( qr(DD)$rank<ncol(DD) & count<1000 ) 
      {
        rand=sample(1:length(TF))
        TF=TF[rand][order(MF[rand])] 
        TB=TB[rand,][order(MF[rand]),]
        DD=crossprod(cbind(TB , NB))
        count=count+1
      }
    if (count==1000) return(TF)   
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
    opt=Optimise(TF,BF,MF,M11,M22,M12,searches,Iterations)
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
    optTF=function(Design,treatlevs,replevs,searches)  {
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
          opt=GenOpt(as.factor(TF),as.factor(Design[,(i+1)]),as.factor(Design[,i]),searches,Iterations[[i]])
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
  
  #**************************** Test for A-optimality *******************************************************
  optEffics=function(TF,BF,treps,breps,ntrts,nblks)   {
    NN= t(table(TF, BF)*(1/sqrt(treps)) ) * (1/sqrt(breps))  
    if (ntrts<=nblks) 
      e=eigen( (diag(ntrts)-crossprod(NN)), symmetric=TRUE, only.values = TRUE)$values[1:(ntrts-1)]     
    else       
      e=c(rep(1,(ntrts-nblks)),eigen((diag(nblks)-tcrossprod(NN)), symmetric=TRUE, only.values = TRUE)$values[1:(nblks-1)])    
    aeff = 1/mean(1/e) 
    deff = exp(sum(log(e))/(ntrts-1))
    list(deff=deff,aeff=aeff)
  }
  
  #******************************************************** A-efficiencies ************************************************************************************
  # A-Efficiencies function 
  A_Efficiencies=function(Design)  {
    strata=ncol(Design)-2
    nunits=nrow(Design)
    treps=tabulate(Design$Treatments)
    ntrts=nlevels(Design$Treatments)
    aeff=rep(1,strata) 
    deff=rep(1,strata)    
    for (i in 1:strata) { 
      nblks=nlevels(Design[,i])
      breps=tabulate(Design[,i])
      if (ntrts>1 & nblks>1) {
        effics=optEffics(Design$Treatments,Design[,i],treps,breps,ntrts,nblks)  
        aeff[i] = effics$aeff       
        deff[i] = effics$deff
      }
    }
    bounds=rep(1,strata)
    blocks=rep(0,strata)
    for (i in 1:strata)   
      blocks[i]=nlevels(Design[,i])    
    if ( all(treps==treps[1]) )
      for (i in 1:strata)  
        if (nunits%%blocks[i] == 0) 
          bounds[i]=upper_bounds(nunits,ntrts,blocks[i])    
    Efficiencies=as.data.frame(cbind(blocks, deff, aeff, bounds))
    colnames(Efficiencies)=c("Blocks","D-Efficiencies","A-Efficiencies", "A-Upper Bounds")
    rnames=c("Main")
    if (strata>1)
      for (i in 1 : (strata-1)) rnames=c(rnames,paste("Sub",i))
    rownames(Efficiencies)=rnames 
    Efficiencies
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
 testInputs=function(treatments,replicates,blocklevels,searches,seed) {  
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
 testout=testInputs(treatments,replicates,blocklevels,searches,seed) 
  if (!isTRUE(testout)) stop(testout)
  if (is.null(seed)) seed=sample(1:100000,1)
  set.seed(seed) 
  # omit any single replicate treatments here 
   treatlevs=treatments[replicates>1]
   replevs = replicates[replicates>1]
  if (is.null(blocklevels)) 
    blocklevels=HCF(replevs)
  nunits=sum(treatlevs*replevs) 
  if (is.null(searches)) 
   searches=min(32, ceiling(4096/nunits))
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
 opt=optTF(Design,treatlevs,replevs,searches) 
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
 
  list(Design=Design,Plan=plan,Incidences=Incidences,Iterations=Iterations,Efficiencies=efficiencies,seed=seed) 
} 