library(SuppDists) #for the rinvGauss (Wald) function

logprofile <- function(z,K,z0,d) log((z-d)/z0)/K

rWALD <- function(n,m,Umeanlog,Usdlog,H,Fmeanlog,Fsdlog,h)
{
  # n: number of random numbers
  # m: wind measurement height
  # Umeanlog: mean log wind speed (m/s) at m
  # Usdlog: sd log wind speed (m/s) at m
  # H: mean plant heigh (m)
  # Fmeanlog: mean log terminal velocity (m/s)
  # Fsdlog: sd log terminal velocity (m/s)
  # h: mean vegetation height (m)
  K <- 0.4    # von Karman constant
  C0 <- 3.125 # Kolmogorov constant
  Aw <- 1.3   # ratio of sigmaw to ustar
  d <- 0.7*h  # zero-plane displacement
  z0 <- 0.1*h # roughness length
  Um <- rlnorm(n,meanlog=Umeanlog,sdlog=Usdlog) # simulate n wind events assuming lognormal distribution
  ustar <- K*Um/log((2-d)/z0)
  U <- (ustar/H)*integrate(logprofile,lower=z0+d,upper=H,K=K,z0=z0,d=d)$value   # compute average wind speed between H and z0+d
  sigma <- sqrt( (4*((Aw)^4)*K*(H-d)/C0) * ustar/U )
  f <- rlnorm(n,meanlog=Fmeanlog,sdlog=Fsdlog) # draw n terminal velocities assuming lognormal distribution
  nu <- H*U/f
  lambda <- (H/sigma)^2
  return(rinvGauss(n,nu=nu,lambda=lambda))
}


distances<-rWALD(n=1000,m=1,Umeanlog=2,Usdlog=.2,H=1.1,Fmeanlog=2,Fsdlog=.5,h=.5)

plot(density(distances), xlim=c(0,10))
