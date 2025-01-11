##     Project Name: Consistent estimation of the proportion via solutions of LS integral equations

### Storey proportion estimator will be applied to discrete p-values
storeyEst <- function(lambda,pvector)  {
   m <- length(pvector)

  over <- m*(1-lambda)
  stunmin <- sum(pvector > lambda)/over 
  st <- min(1,stunmin)

  return(st)
  }
 
### MR estimator 
MREstPiFunc <- function(pvec)
 {
    m = length(pvec)
    pval_sorted = sort(pvec,decreasing = FALSE)

    ## MR estimator with bouding sequence bm = sqrt(2*log(log(m)))/sqrt(m) and dm = sqrt(pval_sorted * (1 - pval_sorted))
    bm = sqrt(2*log(log(m)))/sqrt(m) 
    dm = sqrt(pval_sorted*(1 - pval_sorted))
    pihat = ((1:m)/m - pval_sorted - bm*dm)/(1 - pval_sorted)
    piEst_MR = min(1,max(pihat[2 :(m-2)],0))

    # return other estimators of pi  
    return(piEst_MR)
  }


### Construction I:  Point null: Gaussian
EstPiTypeI <-function(t,x,mu0,str)
{
  n  = length(x)
  xi = seq(-1,1,by=0.01)  # xi are the values to be evaluated by the density function. 0.01 is delta(ksi) for Riemann sum
  w2 = 1/2 + 0*xi      # uniform
  w3 = 1 - abs(xi)     # triangular

  epstmp    = rep(0,length(xi))
  for       (k in 1:length(xi)) {
    epstmp[k] = mean(exp((t*xi[k])^2/2)*cos(t*xi[k]*(x-mu0)))  # change order of summation in the estimator
  }     # above line: heterogeneous nulls and nonnulls          # sum over indices of observations

  if (str=='Uniform') {
          eps0 = sum(w2*epstmp)/100    # 0.01 = 1/100 is the equal space between xi when intergration is
   } else if (  str=='Triangle') {
          eps0 = sum(w3*epstmp)/100
   }   
   return(eps0)
}


############### Construction III ###########################
## the implementation is for Gamma family
EstPiTypeThree = function(t,x,theta0,sigma,Ltrunc)
{
  m  = length(x)
  pm = 0.01
  xi = seq(-1,1,by=pm) 
  L1 = length(xi)
  w3 = 1 - abs(xi)     # use triangular density

  sumMatrix = matrix(0,L1,m)
  epstmp    = rep(0,L1)
  LeftInSeries = 0:Ltrunc 
  for  (k in 1:L1) {
      for (j in 1:m) {
      sumMatrix[k,j] = sum(((-t*xi[k]*x[j])^LeftInSeries)*gamma(sigma)*cos(0.5*pi*LeftInSeries+t*xi[k]/(1-theta0))/(factorial(LeftInSeries)*gamma(sigma+LeftInSeries))) # power series truncated to have summands with indeces 0:Ltrunc
    }
    epstmp[k] = mean(sumMatrix[k,]) # average across sample  
   } 

   ## gather estimate
   eps = 1 - sum(w3*epstmp)*pm
   eps = max(0,eps)
   eps = min(eps,1)
          
   return(eps)
}
 
