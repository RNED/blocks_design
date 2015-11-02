#' @title Efficiency bounds 
#' 
#' @description
#' Finds upper A-efficiency bounds for regular block designs.
#' 
#' @details
#' Upper bounds for the A-efficiency factors of regular nested block designs 
#' (see Chapter 2.8 of John and Williams 1995). Non-trivial bounds
#' are calculated for regular block designs with equal block sizes 
#' and equal replication. All other designs return NA.    
#' 
#' @param nplots the total number of plots in the design.
#' 
#' @param ntrts the total number of treatments in the design.
#' 
#' @param nblocks the total number of blocks in the design.
#' 
#' @references
#' 
#' John, J. A. and Williams, E. R. (1995). Cyclic and Computer Generated Designs. Chapman and Hall, London.
#' 
#' @examples 
#' 
#' # 50 plots, 10 treatments and 10 blocks for a design with 5 replicates and blocks of size 5 
#' upper_bounds(nplots=50,ntrts=10,nblocks=10)
#'
#' @export
upper_bounds=function(nplots,ntrts,nblocks) {
  if (nplots%%ntrts != 0 | nplots%%nblocks != 0 | (ntrts+nblocks-1)>nplots ) return(NA) 
  nreps = nplots/ntrts #replication
  if (nreps%%nblocks == 0) return(1)  
  bsize = nplots/nblocks #block size	
  dual= (ntrts>nblocks & bsize <= ntrts)
  if (dual) {
    temp = nblocks
    nblocks = ntrts
    ntrts = temp
    nreps = nplots%/%ntrts
    bsize = nplots%/%nblocks
  }	
  # average efficiency factor
  ebar =  ntrts*(bsize - 1)/(bsize*(ntrts - 1))
  # average pairwise treatment concurrences within blocks
  lambda = nreps*(bsize - 1)/(ntrts - 1)
  alpha = lambda - floor(lambda)
  bound=ebar
  s2=ntrts*(ntrts-1)*alpha*(1-alpha)/((nreps*bsize)**2) # corrected second moment lower bound

  # this bound is for non-binary designs where bsize>ntrts and can be improved - see John and Williams page 44
  if (bsize > ntrts) {
    k=bsize
    v=ntrts
    b=nblocks
    kp=k%%v
    phi=kp/k
    k0=k%/%v
    U0=1 - kp*(v - kp) / (k*k*(v - 1)) 
    rp=b*kp/v
    lambdap = rp*(kp - 1)/(v - 1)
    alphap = lambdap - floor(lambdap)
    s2p=phi^4*v*(v-1)*alphap*(1-alphap)/((rp*kp)**2) # corrected second moment lower bound
    
    U1= U0 - (ntrts - 2)*s2p/ ( sqrt((ntrts-1)*(ntrts-2)) * (U0*sqrt((ntrts-1)*(ntrts-2)) + (ntrts - 3)*sqrt(s2p)))
    U2= U0 - (1 - U0)*s2p/((1 - U0)*(ntrts - 1) - s2p)
    print(U0)
    print(U1)
    print(U2)
    bound=min(U0,U1,U2,na.rm = TRUE)	
    return(round(bound , 5))		
  }

  if (!isTRUE(all.equal(0,lambda))) {
    if ( alpha < ntrts/(2*(ntrts - 1)) )
      z = alpha*((ntrts + 1)*alpha - 3) else
      z = (1 - alpha)*(ntrts - (ntrts + 1)*alpha)
    s31 = alpha*ntrts*(ntrts - 1)*z/((nreps*bsize)**3)
    if (isTRUE(all.equal(0,alpha))) 
      s32 = alpha*ntrts*(ntrts - 1)*(  (ntrts + 1)*alpha*alpha - 3*alpha - bsize + 2)/((nreps*bsize)**3) else
      s32=s31
    U1= ebar - (ntrts - 2)*s2/ ( sqrt((ntrts-1)*(ntrts-2)) * (ebar*sqrt((ntrts-1)*(ntrts-2)) + (ntrts - 3)*sqrt(s2)))
    U2= ebar - (1 - ebar)*s2/((1 - ebar)*(ntrts - 1) - s2)
    U3= ebar - s2*s2/((ntrts - 1)*(s31+ ebar*s2))	
    U4= ebar - s2*s2/((ntrts - 1)*(s32+ ebar*s2)) 	
    bound=min(U1,U2,U3,U4,na.rm = TRUE)	
  }	
  if (dual) {
    temp = nblocks
    nblocks = ntrts
    ntrts = temp
    bound = (ntrts - 1)/((ntrts - nblocks) + (nblocks - 1)/bound)
  }			
  round(bound,6)
}	
