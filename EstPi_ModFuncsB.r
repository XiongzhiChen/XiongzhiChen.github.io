##     Project Name: Consistent estimation of the proportion via solutions of LS integral equations

#### Bounded Null: Gaussian
EstPi_IV_Gaussian <-function(t,x,a,b)
{
  dxy = 0.01 # norm of partition for Riemman summ
  si = seq(-1,1,by=dxy)  # xi are the values to be evaluated by the density function.
  yi = seq(a,b,by=dxy)   # yi are values for integral wrt y on [a,b]

  lgsi = length(si); lgyi = length(yi)
  aveForSample = matrix(0,nrow=lgyi,ncol=lgsi)  # fixed a yk and sj, sum across sample first
  for  (k in 1:lgyi) {
    for (j in 1:lgsi)  {
     aveForSample[k,j] = (t/(2*pi))*mean(exp((t*si[j])^2/2)*cos(t*si[j]*(x-yi[k])))  # change order of summation in the estimator
   }      # sum over indices of observations
  }
  # sum across yk and sj to get Riemman sum
  eptPi0 = sum(aveForSample)*dxy*dxy
  return(eptPi0)
}

###### one-sided null: Gaussian
EstPi_V_Gaussian <-function(t,x)
{
  dxy = 0.01 # norm of partition for Riemman summ
  si = seq(-1,1,by=dxy)  # xi are the values to be evaluated by the density function.
  yi = seq(0,1,by=dxy)   # yi are values for integral wrt y

  lgsi = length(si); lgyi = length(yi)
  aveForSample = matrix(0,nrow=lgyi,ncol=lgsi)  # fixed a yk and sj, sum across sample first
  for  (k in 1:lgyi) {
    for (j in 1:lgsi)  {
     s1 = t*x*cos(t*si[j]*x*yi[k])*exp((t*si[j]*yi[k])^2/2)
     s2 = t^2*si[j]*yi[k]*sin(t*si[j]*x*yi[k])*exp((t*si[j]*yi[k])^2/2)
     aveForSample[k,j] = (1/(2*pi))*(mean(s1)+mean(s2))  # change order of summation in the estimator
   }      # sum over indices of observations
  }
  # sum across yk and sj to get Riemman sum
  eptPi0 = sum(aveForSample)*dxy*dxy    
  return(eptPi0)
}


## Bounded null: Gamma
EstPi_IV_Gamma = function(t,x,sigma,Ltrunc,a,b)
{
  lgx = length(x)
  dxy = 0.01 # norm of partition for Riemman summ
  si = seq(-1,1,by=dxy)  # xi are the values to be evaluated by the density function.
  yi = seq(a,b,by=dxy)   # yi are values for integral wrt y on [a,b]
  lgsi = length(si); lgyi = length(yi)
  
  LeftInSeries = 0:Ltrunc
  
  aveForSample = matrix(0,nrow=lgyi,ncol=lgsi)  # fixed a yk and sj, sum across sample first
  for  (k in 1:lgyi) {
    for (j in 1:lgsi)  {
    # power series truncated to have summands with indeces 0:Ltrunc 
      sumTmp = double(lgx)
     for (i in 1:lgx) {
      sumTmp[i] = sum(((t*si[j]*x[i]*sigma)^LeftInSeries)*cos(0.5*pi*LeftInSeries-t*si[j]*yi[k])/(factorial(LeftInSeries)*gamma(sigma+LeftInSeries))) 
     }
      aveForSample[k,j] =  mean(sumTmp) # for double integral
   }      # sum over indices of observations
  }
  # sum across yk and sj to get Riemman sum
  eptPi0 = (1/(2*pi))*t*gamma(sigma)*sum(aveForSample)*dxy*dxy    
  return(eptPi0)
}
 

## one-sided null: Gamma
EstPi_V_Gamma = function(t,x,sigma,Ltrunc,b)
{
  lgx = length(x)
  dxy = 0.01 # norm of partition for Riemman summ
  si = seq(-1,1,by=dxy)  # xi are the values to be evaluated by the density function.
  yi = seq(0,1,by=dxy)   # yi are values for integral wrt y
  lgsi = length(si); lgyi = length(yi)
  
  LeftInSeries = 0:Ltrunc
  aveForSample = matrix(0,nrow=lgyi,ncol=lgsi)  # fixed a yk and sj, sum across sample first
  for  (k in 1:lgyi) {
    for (j in 1:lgsi)  {
    # power series truncated to have summands with indeces 0:Ltrunc 
      sumTmp = double(lgx)
     for (i in 1:lgx) {
      sf = sigma*x[i]/gamma(sigma+1+LeftInSeries)-b/gamma(sigma+LeftInSeries)
      sumTmp[i] = sum(((t*si[j]*x[i]*sigma*yi[k])^LeftInSeries)*cos(0.5*pi*LeftInSeries-t*si[j]*yi[k]*b)*sf/factorial(LeftInSeries)) 
     }
      aveForSample[k,j] =  mean(sumTmp) # for double integral
   }      # sum over indices of observations
  }
  # sum across yk and sj to get Riemman sum
  eptPi0 = (1/(2*pi))*t*gamma(sigma)*sum(aveForSample)*dxy*dxy            
  return(eptPi0)
}

