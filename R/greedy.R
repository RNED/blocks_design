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

#********************************************************initialize************************************************************************************
initialize=function(MF,TF,BF) {
  nunits=length(TF)
  mblocks=nlevels(MF)
  nblocks=nlevels(BF)
  repeat {
    Inc=table(BF ,TF)
    N=which.max(Inc)-1 
    i=1+N%%nblocks
    j=1+N%/%nblocks
    main=1+(i-1)%/%(nblocks/mblocks)
    if (Inc[i,j]<=1) break 
    p=sample((1:nunits)[(TF==j)&(BF==i)],1)
    repeat {
      s=sample((1:nunits)[MF==main],1)
      if ( !identical(TF[s],TF[p]) & !identical(BF[s],BF[p])) break
    } 
    TF[c(s,p)]=TF[c(p,s)]
  }
  TF
}

# *****************************************************************

treatlevs=c(7,6,5)
ntrts=sum(treatlevs)
replevs=c(4,3,2)
hcf=HCF(replevs)
TF=NULL
for (i in 1: hcf) 
  TF=c(TF, sample(rep(1:ntrts , rep(replevs/hcf,treatlevs))) )
TF=as.factor(TF)
nunits=sum(treatlevs*replevs)        
blocklevels=12
ntrts=sum(treatlevs)

strata=length(blocklevels)
for (i in 1:strata)
  blocksizes=Sizes(nunits,blocklevels)

facMat= matrix(nrow=prod(blocklevels),ncol=strata)
for (r in 1 : strata) 
  facMat[,r]=gl(prod(blocklevels[1:r]),prod(blocklevels)/prod(blocklevels[1:r])  ) 

Design=facMat[rep(1:length(blocksizes),blocksizes),]
Design=as.data.frame(cbind(rep(1,nunits), Design, rep(1:nunits)))
Design=cbind(Design,as.factor(TF))
Design[]=lapply(Design, factor) 
names(Design)=c("Main","Blocks","Plots","Trts")

table(Design$Blocks ,Design$Trts)
TF=initialize(Design$Main,Design$Trts,Design$Blocks)
table(Design$Blocks ,TF)



