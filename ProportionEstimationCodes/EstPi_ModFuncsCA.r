##     Project Name: Consistent estimation of the proportion via solutions of LS integral equations

# extension of IV for Gaussian
phiFunc = function(x,b,p) {
    if (abs(x) <=b) {
      y1 = (abs(x))^p
     } else {
       y1 = 0
     } 
   return(y1)    
}

#### Bounded Null: Gaussian
EstPi_IVExt_Gaussian <-function(t,x,a,b,p)
{
  dxy = 0.01 # norm of partition for Riemman summ
  si = seq(-1,1,by=dxy)  # xi are the values to be evaluated by the density function.
  yi = seq(a,b,by=dxy)   # yi are values for integral wrt y on [a,b]

  lgsi = length(si); lgyi = length(yi)
  aveForSample = matrix(0,nrow=lgyi,ncol=lgsi)  # fixed a yk and sj, sum across sample first
  for  (k in 1:lgyi) {
    for (j in 1:lgsi)  {
     aveForSample[k,j] = (t/(2*pi))*mean(phiFunc(yi[k],b,p)*exp((t*si[j])^2/2)*cos(t*si[j]*(x-yi[k])))  # change order of summation in the estimator
   }      # sum over indices of observations
  }
  # sum across yk and sj to get Riemman sum
  eptPi0 = sum(aveForSample)*dxy*dxy
  return(eptPi0)
}

## sim data
GenGaussianExtA_IV = function(betaMin,betaMax,a,b,p,minGap) {
    # beta in (a,b)
    betavecNull = runif(m0,a,b)  
    betavecAltA = runif(floor(0.5*m1),b+minGap,b+betaMax)
    betavecAltB = runif(m1-floor(0.5*m1),a+betaMin,a-minGap)
    betavec = c(c(betavecNull,betavecAltA),betavecAltB)
    # compute new prop
    Tmp = double(m0)
    for (ib in 1:m0) {
       Tmp[ib]= phiFunc(betavecNull[ib],b,p)
     }
     phiPro = sum(Tmp)/(m0+m1)   
    return(list(betavec,phiPro))
  } 