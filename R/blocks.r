#' @title Block designs 
#' 
#' @description
#' 
#' Constructs nested block designs for unstructured treatment sets with arbitrary replication, not necessarily all equal, and arbitrary depth of nesting.
#' 
#' @details
#' 
#' The \code{blocks} function constructs optimized nested block designs for unstructured treatment sets where treatments can have any arbitrary replication, not necessarily all equal, 
#' and blocks can have any feasible depth of nesting.
#' 
#' The \code{treatments} list contains the total number of treatments in the design partitioned into
#' sets of equally replicated treatments where the replication can vary between sets but must be constant within sets.  
#' The \code{replicates} list provides the replication level for each treatment set and must be the same length as 
#' the \code{treatments} list. The individual treatments are labelled consecutively according to the ordering of the sets
#' in the two lists and different treatment sets with the same replication can be used, if required.
#'  
#' The \code{blocklevels} list defines the blocks structure of the design where the first level in the list is the number of main blocks 
#' and the succesive levels, if any, are the levels of sub-blocks in a hierarchy of nested sub-blocks. Each nested level is the number of sub-blocks
#' nested in each block of the preceding stratum. The length of the list is the number of strata
#' in the design and the running product of the first k-levels of the \code{blocklevels} list is the total numbers of 
#' blocks in the kth stratum. The average block size in any stratum is the total number of blocks in that stratum divided by 
#' the total number of plots therefore the block sizes are exactly equal if the quotient is an integer. 
#' Otherwise, the block sizes differ by, at most, a single unit. 
#' 
#' The default design is a main blocks design with a complete set of othogonal blocks. 
#' 
#' General block designs are constructed algorithmically by a swapping algorithm that seeks to maximize the determinant of
#' the information matrix (D-optimality). Beginning with the main blocks stratum, pairs of treatments are swapped between
#' sub-blocks nested within blocks until a local optima is attained. If the number of searches
#' is greater than one, the algorithm makes a number of random swaps and then continues with the improving swaps
#' until another local optima is attained. After the required number of searches, the best overall design
#' for that stratum is restored and the process is repeated for the next stratum in the hierarchy. Eventually the bottom stratum
#' is optimized after which the algorithm stops.  
#' 
#' Special lattice block designs for v**2 equally replicated treatments in blocks of size v with k replicates 
#' are constructed algebraicaly when k <= 3 or when v is prime or prime-power and k <= v+1 or when v = 10 and k <= 4. 
#' The \code{crossdes} package is required for lattice designs with prime-power v.
#'   
#' Optimized designs are fully randomized with each set of nested blocks randomized within the preceding set of blocks and with
#' the  treatments randomized within the bottom set of blocks.
#'  
#' @param treatments a list partitioning the total number of treatments into equally replicated treatment sets.   
#' 
#' @param replicates a list assigning a replication level to each equally replicated treatment set. 
#' 
#' @param blocklevels a list of nested block levels where the first level is the number of main blocks
#' and the remaining levels, if any, are the levels of sub-blocks in  a hierarchy of nested sub-blocks.
#' The default setting is the highest common factor of the \code{replicates} list.
#'  
#' @param searches the number of local optima searched during optimization. The default setting is the minimum of 64 or the integer quotient of 4096 
#' divided by the number of plots.
#' 
#' @param seed an integer seed for initializing the random number generator where a design must be reproducible. The default 
#' setting is a random seed.
#' 
#' @return  
#' \item{Design}{Data frame showing the listing of treatments allocated to blocks}
#' \item{Plan}{Data frame showing a plan of treatments allocated to sub-plots within blocks}
#' \item{Incidences}{List of blocks-by-treatments incidence matrices, one for each stratum in the design}
#' \item{Efficiencies}{Data frame showing the A-efficiency factor for each stratum in the design together with an upper bound, where available}
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
#' @export
#' 

blocks = function(treatments, replicates, blocklevels=NULL, searches=NULL, seed=NULL) { 
  
  #********************************************************Block Sizes ************************************************************************************
  Sizes=function(mainSizes,slevs) {
    if (max(mainSizes)==min(mainSizes) ) {
      bsize=mainSizes[1] %/% slevs
      resid=mainSizes[1] %% slevs
      newsizes=rep( rep(bsize, slevs), length(mainSizes) )  +  rep(  c( rep(1, resid), rep(0, (slevs-resid)  ) ),   length(mainSizes)         )   
    } else {
      newsizes=vector(length=slevs*length(mainSizes))
      for (z in 1: length(mainSizes)) {
        bsize = mainSizes[z] %/% slevs
        resid = mainSizes[z] %% slevs
        for (i in 1 : slevs) {
          newsizes[(z-1)*slevs+i] = bsize + (resid>0)
          resid=resid-1
        }
      }   
    }
    newsizes 
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
      BM[(as.numeric(MF)==i),]=scale(BM[(as.numeric(MF)==i),], center = TRUE, scale = FALSE) # scaling within main blocks
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
  
  #********************************************************Determinants of swaps ************************************************************************************
  DetMat=function(M11,M22,M12,Samp,TF,RF) {
    detmat=vector("list", length(Samp))
    for (i in 1:length(Samp)) {
      trts=as.numeric(TF[Samp[[i]]])
      rblks=as.numeric(RF[Samp[[i]]])
      TT=M11[trts,trts,drop=FALSE]
      BB=M22[rblks,rblks,drop=FALSE]
      TB=M12[trts,rblks,drop=FALSE]
      TT=TT-crossprod(t(diag(TT)),t(rep(1,ncol(TT))))
      TT=TT+t(TT)
      BB=BB-crossprod(t(diag(BB)),t(rep(1,ncol(TT))))
      BB=BB+t(BB)
      TB=TB-crossprod(t(diag(TB)),t(rep(1,ncol(TT))))
      TB=1+TB+t(TB)
      detmat[[i]]=TB**2-TT*BB
    }
    detmat
  }
  
  #******************************************************** Updates variance matrix ************************************************************************************
  UpDate=function(M11,M22,M12,si,sj,TF,RF) {
    M11T = M11[TF[si],]-M11[TF[sj],]
    M12T = M12[TF[si],]-M12[TF[sj],]
    M12B = M12[,RF[si]]-M12[,RF[sj]]
    M22B = M22[,RF[si]]-M22[,RF[sj]]
    # updating vectors
    m11=M11T[TF[si]]-M11T[TF[sj]]
    m22=M22B[RF[si]]-M22B[RF[sj]]
    m12=M12T[RF[si]]-M12T[RF[sj]]
    f = sqrt(2+m11+m22-2*m12)
    m = f/sqrt(1-2*m12-m11*m22+m12*m12)/2
    Z1 = (M12B-M11T)/f 
    Z2 = (M22B-M12T)/f 
    W1 = (M11T+M12B - Z1*(m22-m11)/f)*m
    W2 = (M12T+M22B - Z2*(m22-m11)/f)*m
    M11 = M11 - crossprod(t(Z1)) + crossprod(t(W1))
    M22 = M22 - crossprod(t(Z2)) + crossprod(t(W2))
    M12 = M12 - crossprod(t(Z1),t(Z2)) + crossprod(t(W1),t(W2))
    up=list(M11=M11,M22=M22,M12=M12)
  } # end of function
 
  #******************************************************** General optimization ************************************************************************************
  Optimise=function(TF,BF,MF,M11,M22,M12,searches)   {
    # first stage finds an optima by optimizing samples of increasing size in powers of 2
    globrelD=1
    locrelD=1
    globTF=TF
    nunits=length(TF)
    mainSets=split(rep(1:nunits),MF)
    Samp=vector("list", nlevels(MF)) 
    for (r in 1 : searches) {
      nSamp=ceiling(  min(nunits,36)*tabulate(MF)/nunits) #sample size is smallest of 36 or nunits or has at least one sample per restriction
      repeat {
        relD=0
        for (i in 1:nlevels(MF)) Samp[[i]]=sort(sample(mainSets[[i]],nSamp[i]))						
        detmat=DetMat(M11,M22,M12,Samp,TF,BF)
        grelD=1
        for (i in 1:nlevels(MF)) {
          N=which.max(detmat[[i]])
          ti=1+(N-1)%%nrow(detmat[[i]])
          tj=1+(N-1)%/%nrow(detmat[[i]])
          relD=detmat[[i]][ti,tj]
          if (relD>grelD) {
            gsi=Samp[[i]][ti]
            gsj=Samp[[i]][tj]
            grelD=relD
          }	
        }			
        if (grelD>1.00000001) {
          up=UpDate(M11,M22,M12,gsi,gsj,TF,BF)
          M11=up$M11
          M22=up$M22
          M12=up$M12
          TF[c(gsi,gsj)]=TF[c(gsj,gsi)]
          locrelD=grelD*locrelD				
        } else if ( sum(nSamp) < min(nunits,512)) {
          nSamp=2*nSamp # doubles sample size
          if (sum(nSamp)>nunits) nSamp=tabulate(MF) # ensures sample size not greater than the population
        } else break
      } # repeat	
      if (locrelD>globrelD) {
        globrelD=locrelD
        globTF=TF
      }
      if (searches>r) {
        # escape local optima
        prop_change=1
        for (iswap in 1 : 6) {
          icount=0
          dswap=1
          while (icount<100 & (dswap<0.01 | dswap>.999)) {
            icount=icount+1
            s=sample(rep(1:nunits)[MF==sample(nlevels(MF),1)],2)
            # calculates the proportional change in the determinant of the design information due to swapping treatments on plots s1 and s2
            if ( (TF[s[1]]!=TF[s[2]]) & (BF[s[1]]!=BF[s[2]]) ) 	
              dswap=(1+M12[TF[s[1]],BF[s[2]]]+M12[TF[s[2]],BF[s[1]]]-M12[TF[s[1]],BF[s[1]]]-M12[TF[s[2]],BF[s[2]]])**2-
              (2*M11[TF[s[1]],TF[s[2]]]-M11[TF[s[1]],TF[s[1]]]-M11[TF[s[2]],TF[s[2]]])*
              (2*M22[BF[s[1]],BF[s[2]]]-M22[BF[s[1]],BF[s[1]]]-M22[BF[s[2]],BF[s[2]]])
          }
          #updates matrices
          if (icount<100) {
            prop_change=prop_change*dswap
            up=UpDate(M11,M22,M12,s[1],s[2],TF,BF)  
            M11=up$M11
            M22=up$M22
            M12=up$M12
            TF[c(s[1],s[2])]=TF[c(s[2],s[1])]	
          }
        } 
        locrelD=locrelD*prop_change
      } 
    } # next for
    globTF
  } # end of function
  
  #******************************************************** fullrank randomised design  ************************************************************************************
  FullRank=function(MF,TF,TB,NB,DD) {
    count=0
    while ( (qr(DD)$rank<ncol(DD))&(count<1000) ) {
      rand=sample(1:length(TF))
      TF=TF[rand][order(MF[rand])]
      TB=TB[rand,][order(MF[rand]),]
      DD=crossprod(cbind(TB,NB))
      count=count+1
    } # end of while
    fullrank=list(TF=TF,DD=DD,count=count)	
  }
  
  #******************************************************** Initializes design***************************************************************************************
  GenOpt=function(TF,NF,MF,searches)  {
    singular=FALSE
    T=TreatContrasts(MF,TF)
    N=BlockContrasts(MF,NF) 
    DD=crossprod(cbind(T,N))
    count=0
    # attempts to find a full rank randomisation	
    if (qr(DD)$rank<ncol(DD)) {
      fullrank=FullRank(MF,TF,T,N,DD)
      TF=fullrank$TF
      count=fullrank$count
      DD=fullrank$DD
    }
    if (count<100) {
      V=chol2inv(chol(DD))
      M11=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))	
      M22=matrix(0,nrow=nlevels(NF),ncol=nlevels(NF))
      M12=matrix(0,nrow=nlevels(TF),ncol=nlevels(NF))
      M11[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[1:(nlevels(TF)-1),1:(nlevels(TF)-1),drop=FALSE]
      M12[1:(nlevels(TF)-1),1:(ncol(V)-nlevels(TF)+1)]=V[1:(nlevels(TF)-1),nlevels(TF):ncol(V),drop=FALSE]
      M22[1:(ncol(V)-nlevels(TF)+1),1:(ncol(V)-nlevels(TF)+1)]=V[nlevels(TF):ncol(V),nlevels(TF):ncol(V),drop=FALSE]
      sortR=c(1:nlevels(NF))%%(nlevels(NF)/nlevels(MF))==0
      sortC=NULL
      sortN=NULL
      perm=order(order(c(sortR,sortC,sortN)))
      M12=M12[,perm]
      M22=M22[perm,perm]	
      TF=Optimise(TF,NF,MF,M11,M22,M12,searches)
    }	
    TF
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
 
  #******************************************************** builds lattice designs  ************************************************************************************
  optTF=function(Design,treatlevs,replevs,searches) {
    nunits=nrow(Design)
    ntrts=sum(treatlevs)
    hcf=HCF(replevs)
    v=sqrt(ntrts)
    strata=ncol(Design)-1
    orthbsize=nunits/hcf 
    ortho=0
    for (i in 1 : strata) 
      if (all( tabulate(Design[,i+1]) %% orthbsize == 0)) ortho=i else  break
    reglat=(max(replevs)==min(replevs) & max(Design[,i+1])==v*replevs[1] & identical(v,ntrts%/%v) )
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
          TF=GenOpt(as.factor(TF),as.factor(Design[,(i+1)]),as.factor(Design[,i]),searches)
        }
      }
    }
    TF   
  }
  
  #******************************************************** Randomizes blocks************************************************************************************
  randBlocks=function(Design) {
    #randomise blocks within strata
    Design=as.data.frame( cbind( Design[,1:(ncol(Design)-1)], rep(1:nrow(Design)), Design[,(ncol(Design))] ) )   
    Design[]=lapply(Design, factor)  
    for (r in 2 : (ncol(Design)-1)){
      levels( Design[,r])=sample(nlevels( Design[,r]) )
      Design[,r]=as.numeric(levels(Design[,r]))[Design[,r]]
    }
    # re-order all numeric columns
    Design=Design[ do.call(order, Design), ]
    Design=Design[,c(1:(ncol(Design)-2),ncol(Design))]
    Design
  }
 
  #******************************************************** Puts back single rep treatments *************************************************************************** 
  fullDesign=function(Design,treatments,replicates,blocksizes,blocklevels) {
    strata=ncol(Design)-2
    TF=as.factor(Design[,ncol(Design)])
    nunits=sum(treatments*replicates)  
    addTF=c( (nlevels(TF)+1 ) : sum(treatments))
    TF=as.factor(c(TF,addTF[sample(length(addTF))])) 
    trtlabs=NULL  
    extlabs=NULL
    index=0
    for (i in 1 : length(treatments)) {
      if (replicates[i]>1) 
        trtlabs=c(trtlabs, rep( (index+1):(index+treatments[i]) ))
      else extlabs=c(extlabs, rep( (index+1):(index+treatments[i]) ))
        index=index+treatments[i] 
    }
    trtlabs=c(trtlabs,extlabs)
    levels(TF)=trtlabs
    TF=as.numeric(levels(TF))[TF]
    BF=c( rep( 1:length(blocksizes),blocksizes))
    newblocksizes=nunits
    for (i in 1:strata)
     newblocksizes=Sizes(newblocksizes,blocklevels[i])   
    BF=c(BF,  rep( 1:length(blocksizes),(newblocksizes-blocksizes) ) )
    # full TF in blocks
    TF=TF[order(BF)]
    Design = matrix(1,nrow=nunits,ncol=(strata+1))
    for (r in 1 : strata) 
      Design[,r+1]=rep(facMat[,r],newblocksizes)
    Design=cbind(Design,as.factor(TF))  
    Design
  }  
  
  #******************************************************** A-efficiencies ************************************************************************************
  # A-Efficiencies function 
  A_Efficiencies=function(Design)  {
    strata=ncol(Design)-1
    nunits=nrow(Design)
    TF=Design$Treatments
    treps=tabulate(TF)
    ntrts=nlevels(TF)
    aeff=rep(0,strata)  
    for (i in 1:strata) { 
      bSize=tabulate(Design[,i])
      nblks=nlevels(Design[,i])
      if (ntrts<=nblks) {
        X=crossprod( diag(1/sqrt(bSize),nrow = nblks),  table(Design[,i],TF ) )
        A= diag(ntrts) - crossprod(crossprod(t(X), diag(1/sqrt(treps),nrow = ntrts)))   
        aeff[i] = 1/mean(1/eigen(A, symmetric=TRUE, only.values = TRUE)$values[1:(ntrts-1)])
      } else {
        X=crossprod( diag(1/sqrt(treps),nrow = ntrts),  table(TF,Design[,i] ) )
        A=diag(nblks) - crossprod(crossprod(t(X), diag(1/sqrt(bSize), nrow = nblks)))  
        aeff[i]=1/mean(1/eigen(A, symmetric=TRUE, only.values = TRUE)$values[1:(nblks-1)])    
        aeff[i]=(ntrts-1)/ (ntrts-nblks+(nblks-1)/aeff[i] )
      }
    }
    bounds=rep(NA,strata)
    blocks=rep(0,strata)
    if (max(treps)==min(treps))
      for (i in 1:strata)   
        bounds[i]=upper_bounds(nunits,ntrts,nlevels(Design[,i]))  
    for (i in 1:strata)   
      blocks[i]=nlevels(Design[,i])
    Efficiencies=as.data.frame(cbind(blocks, aeff, bounds))
    colnames(Efficiencies)=c("Blocks","A-Efficiencies", "Upper Bounds")
    rnames=c("Main")
    if (strata>1)
      for (i in 1 : (strata-1)) rnames=c(rnames,paste("Sub",i))
    rownames(Efficiencies)=rnames 
    Efficiencies
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
    return(" Treatments and replicates must contain only finite integers ")
  if ( length(treatments)!=length(replicates) ) 
    return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))
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
  if (  sum(treatments*replicates) < (prod(blocklevels) + sum(treatments)-1) ) 
    return("Design cannot be fitted :  too many blocks and treatments for the available plots")  
  return(TRUE)
  }
 
 #******************************************************** Plan output************************************************************************************
 Plan=function(Design,facMat)  {
   strata=ncol(Design)-1
   bSizes=tabulate(Design[,strata])
   plotTrts=matrix(nrow=length(bSizes),ncol=max(bSizes)) 
   counter=0
   for (i in 1:length(bSizes)) {
       plotTrts[i,c(1 : bSizes[i])]=Design[ c((1+counter) : (counter+bSizes[i]))  , strata+1]  
       counter=counter+bSizes[i]
   }
   plotTrts[is.na(plotTrts)]  = " "
   Plan=as.data.frame(cbind(facMat,rep(" ",length(bSizes)),plotTrts))
   designnames=c(designnames,"Sub_plots")
   for (i in 1:ncol(plotTrts))
     designnames=c(designnames,i)
   colnames(Plan)=designnames
   Plan
 }
 
 #********************************************************builds design ************************************************************************************ 
 testout=testInputs(treatments,replicates,blocklevels,searches,seed) 
  if (testout!=TRUE) return(testout)
  if (is.null(seed)) seed=sample(1:100000,1)
  set.seed(seed) 
  # omit any single replicate treatments here 
   if (all(replicates==1)) {
    treatlevs=treatments
    replevs = replicates
  } else {
   treatlevs=treatments[replicates>1]
   replevs = replicates[replicates>1]
  }
  if (is.null(blocklevels)) 
    blocklevels=HCF(replevs)
  nunits=sum(treatlevs*replevs) 
  if (is.null(searches)) 
   searches=min(64, floor(4096/nunits))
  if (all(blocklevels==1))
    blocklevels=1
  else
  blocklevels=blocklevels[blocklevels>1]
  strata=length(blocklevels)	
  cumblocklevs=1
  for (i in 1 : strata) 
    cumblocklevs=c(cumblocklevs,cumblocklevs[i]*blocklevels[i])
  blocksizes=nunits
  for (i in 1:strata)
    blocksizes=Sizes(blocksizes,blocklevels[i])
  facMat= matrix(nrow=cumblocklevs[strata+1],ncol=strata)
  for (r in 1 : strata) 
    facMat[,r]=gl(cumblocklevs[r+1],cumblocklevs[strata+1]/cumblocklevs[r+1])  
 Design= matrix(1,nrow=nunits,ncol=(strata+1))
  for (r in 1 : strata) 
    Design[,r+1]=rep(facMat[,r],blocksizes)
  TF=optTF(Design,treatlevs,replevs,searches)
  Design=cbind(Design,as.factor(TF))   
  # add back single replicate treatments here 
  if ( !all(replicates>1) & !all(replicates==1) ) 
   Design= fullDesign(Design,treatments,replicates,blocksizes,blocklevels)  
  Design=randBlocks(Design)
  # arrange blocksizes in new order
  blocks=as.factor(Design[,(ncol(Design)-1)])
  levels(blocks)=order(unique(blocks))
  blocksizes=tabulate(as.numeric(levels(blocks))[blocks])
  # design matrix for new TF where block sizes are not necessarily all equal 
  for (r in 1 : strata) 
    Design[,(r+1)]=rep(facMat[,r],blocksizes)
  Design[]=lapply(Design, factor)	
  # drop leading column
  Design=Design[,c(2:ncol(Design))]
  designnames=c("Main_blocks")
  if (strata>1)
    for (i in 1:(strata-1))
      designnames=c(designnames,paste("Sub",i,"_blocks", sep=""))
  colnames(Design)=c(designnames,"Treatments")   
  rownames(Design) = NULL 
  # Incidence matrix for each stratum
  Incidences=vector(mode = "list", length =strata )
  for (i in 1:strata)
    Incidences[[i]]=table( Design[,i] ,Design[,strata+1]) 
 
  list(Design=Design,Plan=Plan(Design,facMat),Incidences=Incidences,Efficiencies=A_Efficiencies(Design),seed=seed)
} 