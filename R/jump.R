# Jump-Diffusion Models
#  JumpDiffusionMerton('p',110,109,.25,0.05,.2,5,.6)  # 3.1219
#  greeks.JDB(JumpDiffusionBates('p',100,100,0.5,0.05,0.05,0.2,5,-0.05,.1))

JumpDiffusionMerton <- function(type=c("call","put"),S,X,T,r,v,lambda,gamma) {
  delta <- sqrt(gamma * v^2 / lambda)
  Z <- sqrt(v^2 - lambda * delta^2)
  s <- 0
  i <- 0:50
  vi <- sqrt(Z^2 + delta^2 * (i/T))
  if(match.arg(type) == "call")
    value <- sum(exp(-lambda*T)*(lambda*T)^i / factorial(i) * call.value(S,X,r,r,T,vi)$value)
  else
    value <- sum(exp(-lambda*T)*(lambda*T)^i / factorial(i) * put.value(S,X,r,r,T,vi)$value)
  structure(list(value=value,type=type,S=S,X=X,T=T,r=r,v=v,lambda=lambda,gamma=gamma), class="JDM")
}

greeks.JDM <- function(x, dS=0.01, ...) {
  JDMplus <- with(x, JumpDiffusionMerton(type,S+dS,X,T,r,v,lambda,gamma))$value
  JDMminus <- with(x, JumpDiffusionMerton(type,S-dS,X,T,r,v,lambda,gamma))$value
  delta <- (JDMplus - JDMminus) / (2*dS)
  elasticity <- (JDMplus - JDMminus) / (2*dS) * x$S / x$value
  gamma <- (JDMplus - 2*x$value + JDMminus) / dS^2
  gammaP <- x$S / 100*gamma

  JDMV0.01plus <- with(x, JumpDiffusionMerton(type,S,X,T,r,v+0.01,lambda,gamma))$value
  JDMV0.01minus <- with(x, JumpDiffusionMerton(type,S,X,T,r,v-0.01,lambda,gamma))$value
  vega <- (JDMV0.01plus - JDMV0.01minus) / 2
  vegaP <- x$v / 0.1 * vega

  rho <- (with(x, JumpDiffusionMerton(type,S,X,T,r+0.01,v,lambda,gamma))$value - 
          with(x, JumpDiffusionMerton(type,S,X,T,r-0.01,v,lambda,gamma))$value ) / 2
  list(Delta=delta,Elasticity=elasticity, Gamma=gamma, GammaP=gammaP, Vega=vega, VegaP=vegaP, Rho=rho)
}

JumpDiffusionBates <- function(type=c("call","put"),S,X,T,r,b,v,lambda,avgk,delta) {
  gam0 <- log(1+avgk)
  i <- 0:50
  bi <- b - lambda*avgk + gam0*(i/T)
  vi <- sqrt(v^2 + delta^2 * (i/T))
  if(match.arg(type)=="call")
    value <- sum(exp(-lambda*T) * (lambda*T)^i / factorial(i) * call.value(S,X,bi,r,T,vi)$value)
  else
    value <- sum(exp(-lambda*T) * (lambda*T)^i / factorial(i) * put.value(S,X,bi,r,T,vi)$value)
  structure(list(value=value,type=type,S=S,X=X,T=T,r=r,b=b,v=v,lambda=lambda,avgk=avgk,delta=delta), class="JDB")
}

greeks.JDB <- function(x, dS=0.01, ...) {
  JDBplus <- with(x, JumpDiffusionBates(type,S+dS,X,T,r,b,v,lambda,avgk,delta))$value
  JDBminus <- with(x, JumpDiffusionBates(type,S-dS,X,T,r,b,v,lambda,avgk,delta))$value
  delta <- (JDBplus - JDBminus) / (2*dS)
  elasticity <- (JDBplus - JDBminus) / (2*dS) * x$S / x$value
  gamma <- (JDBplus - 2*x$value + JDBminus) / dS^2
  gammaP <- x$S / 100*gamma

  JDBV0.01plus <- with(x, JumpDiffusionBates(type,S,X,T,r,b,v+0.01,lambda,avgk,delta))$value
  JDBV0.01minus <- with(x, JumpDiffusionBates(type,S,X,T,r,b,v-0.01,lambda,avgk,delta))$value
  vega <- (JDBV0.01plus - JDBV0.01minus) / 2
  vegaP <- x$v / 0.1 * vega

  rho <- (with(x, JumpDiffusionBates(type,S,X,T,r+0.01,b,v,lambda,avgk,delta))$value - 
          with(x, JumpDiffusionBates(type,S,X,T,r-0.01,b,v,lambda,avgk,delta))$value ) / 2
  list(Delta=delta,Elasticity=elasticity, Gamma=gamma, GammaP=gammaP, Vega=vega, VegaP=vegaP, Rho=rho)
}
