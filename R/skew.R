# Skew and Kurtosis Adjustments

JarrowRuddSkewKurtosis <- function(type=c("call","put"),
                                   S, #stock=42,
                                   X, #strike=40,
                                   r, #r=0.1,
                                   b, #b=0.1,
                                   v, #sigma=0.2,
                                   Time, #T=0.5,
                                   skew=0,
                                   kurt=3) 
{
  stock <- S
  strike <- X
  sigma <- v
  T <- Time

  d1 <- ( log(stock/strike) + (r + sigma^2/2)*T )/ (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  aX <- (strike * sigma * sqrt(T*2*pi))^-1 * exp(-d2^2 / 2)
  daX <- aX * (d2 - sigma * sqrt(T)) / (strike* sigma * sqrt(T))
  daXX <- aX / (strike^2 * sigma * sqrt(T)) * ((d2 - sigma * sqrt(T))^2 - sigma * sqrt(T) * (d2 - sigma * sqrt(T)) - 1)
  q <- sqrt(exp(sigma^2 * T) - 1)
  GA <- 3*q*q^3
  gAA <- 16 * q^2 + 15 * q^4 + 6 * q^6 + q^8 + 3
  Lambda1 <- skew - GA
  Lambda2 <- kurt - gAA

  Q3 <- -(stock * exp(r * T))^3 * (exp(sigma^2 * T) - 1)^(3/2) * exp(-r*T) / 6 * daX
  Q4 <- (stock * exp(r*T))^4 * (exp(sigma^2*T) - 1) ^ 2 * exp(-r*T) / 24 * daXX

  cp <- BlackSholesMerton('c',stock,strike,b,r,T,sigma)$value + Lambda1*Q3 + Lambda2*Q4
  if(match.arg(type) == "call")
    cp
  else
    cp - stock*exp((b-r)*T) + strike*exp(-r*T)
}

CorradoSuSkewKurtosis <- function(type=c("call","put"),
                                   S, #stock=42,
                                   X, #strike=40,
                                   r, #r=0.1,
                                   b, #b=0.1,
                                   v, #sigma=0.2,
                                   Time, #T=0.5,
                                   skew=0,
                                   kurt=3) 
{
  stock <- S
  strike <- X
  sigma <- v
  T <- Time

  d1 <- ( log(stock/strike) + (b + sigma^2/2)*T )/ (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  Q4 <- 1/24 * stock * sigma * sqrt(T) * ((d1^2-1-3*sigma*sqrt(T)*d2) * dnorm(d1) + sigma^3 * T^1.5 * pnorm(d1))
  Q3 <- 1/6 * stock * sigma * sqrt(T) * ((2*sigma*sqrt(T)-d1) * dnorm(d1) + sigma^2 * T * pnorm(d1))
  
  cp <- call.value(stock,strike,b,r,T,sigma)$value + skew*Q3 + (kurt-3)*Q4
  if(match.arg(type) == "call")
    cp
  else
    cp - stock*exp((b-r)*T) + strike*exp(-r*T)  # put via put-call parity
}

ModifiedCorradoSuSkewKurtosis <- function(type=c("call","put"),
                                   S, #stock=42,
                                   X, #strike=40,
                                   r, #r=0.1,
                                   b, #b=0.1,
                                   v, #sigma=0.2,
                                   Time, #T=0.5,
                                   skew=0,
                                   kurt=3) 
{
  stock <- S
  strike <- X
  sigma <- v
  T <- Time

  w <- skew/6 * sigma^3 * T^1.5 + kurt/24 * sigma^4 * T^2
  d <- (log(stock/strike) + (b+sigma^2/2)*T - log(1+w)) / (sigma * sqrt(T))
  Q3 <- 1/(6*(1+w))*stock*sigma*sqrt(T)*(2*sigma*sqrt(T)-d)*dnorm(d)
  Q4 <- 1/(24*1+w)*stock*sigma*sqrt(T)*(d^2-3*d*sigma*sqrt(T)+3*sigma^2*T-1)*dnorm(d)

  cp <- call.value(stock,strike,b,r,T,sigma)$value + skew*Q3 + (kurt-3)*Q4
  if(match.arg(type) == "call")
    cp
  else
    cp - stock*exp((b-r)*T) + strike*exp(-r*T)  # put via put-call parity
}


RubinsteinSKB <- function(type,S,X,T,r,b,v,skew,kurt,n,expansion=c('edgeworth','gram-charlier')) {
  z <- ifelse(grepl("c",match.arg(type)), 1, -1)
  dt <- T/n
  u <- exp( (b-v^2 / 2) * dt + v*sqrt(dt))
  d <- exp( (b-v^2 / 2) * dt - v*sqrt(dt))
}
