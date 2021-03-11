#### TODO
### helpers to move:



fischer.skewness <- function(r) {
  # r: lognormal asset returns
  
}  

CRRtree <- function(type=c('ac','ap','ec','ep'),
                    S, #stock=42,
                    X, #strike=40,
                    r, #r=0.1,
                    b, #b=0.1,
                    v, #sigma=0.2,
                    Time, #T=0.5,
                    N=52, lambda=1, drift=0,...) 
{
  # e.g. from Haug
  # S=42, X=40, r=0.1, b=0.1, v=0.2, Time=0.5, N=52
  # http://www.sitmo.com/article/binomial-and-trinomial-trees/
  stock <- S
  strike <- X
  sigma <- v
  T <- Time

  type <- match.arg(type)
  z <- ifelse(type %in% c('ac','ec'), 1, -1)
  american <- ifelse(type %in% c('ac','ap'), TRUE, FALSE)

  # Example from Haug PlainTrees.xls
  dt <- T/N

  if(lambda != 1 && drift != 0)
    warning(paste('only one of',sQuote('lambda'),'or',sQuote('drift'),'should be set'))

  #  Generalized CRR allows for a lambda 'tilt' parameter
  #  http://www.goddardconsulting.ca/option-pricing-binomial-alts.html#crrdrift
  #  see also Chung and Shih (2007) for G-CRR extension
  if(lambda != 1)
    warning(paste(sQuote('lambda'),'parameter is experimental.'))

  u <- exp( (drift*dt) + lambda * sigma * sqrt(dt))
  #d <- 1/u
  d <- exp( (drift*dt) -(1/lambda)*sigma*sqrt(dt))
  p = (exp(b * dt) - d) / (u - d)
  dis <- exp(-r*dt)

  # stock price
  S <- matrix(0,ncol=N+1,nrow=N+1) 
  S <- lapply(0:N, function(N.) stock * u^(0:N.) * d^(N.:0))
 
  # option price
  P <- vector("list", length(S))
  P[[length(P)]] <- ifelse(z* (S[[length(S)]]-strike) > 0, z * (S[[length(S)]]-strike), 0)
  if(american) {
    for(i in N:1)
      P[[i]] <- mapply(max, dis*(p*P[[i+1]][1+(1:i)] + (1-p)*P[[i+1]][1:i]),z*(S[[i]]-strike))
  } else {
    for(i in N:1)
      P[[i]] <- dis*(p*P[[i+1]][1+(1:i)] + (1-p)*P[[i+1]][1:i])
  }
  g <- list()
  g$delta <- diff(P[[2]]) / (S[[1]]*u - S[[1]]*d)
  g$gamma <- (diff(P[[3]][2:3]) / (stock*u^2-stock) - diff(P[[3]][1:2]) / (stock-stock*d^2)) / (0.5*(stock*u^2-stock*d^2))
  g$theta <- (P[[3]][2] - P[[1]]) / (2*dt) / 365
  structure(list(Stock=S,Price=P,Greeks=g), class=c("binomialtree", "CRRtree"))
}

JRtree <- function(type=c('ac','ap','ec','ep'),
                   S, #stock=42,
                   X, #strike=40,
                   r, #r=0.1,
                   b, #b=0.1,
                   v, #sigma=0.2,
                   Time, #T=0.5,
                   N=52, ...) 
{
  stock <- S
  strike <- X
  sigma <- v
  T <- Time

  type <- match.arg(type)
  z <- ifelse(type %in% c('ac','ec'), 1, -1)
  american <- ifelse(type %in% c('ac','ap'), TRUE, FALSE)
  dt <- T/N
  #u <- exp(sigma * sqrt(dt))
  #d <- exp(-sigma * sqrt(dt))
  u <- exp( (b-0.5*sigma^2)*dt + sigma*sqrt(dt))
  d <- exp( (b-0.5*sigma^2)*dt - sigma*sqrt(dt))

  dis <- exp(-r*dt)
  p <- 0.5 #*(1.0 + (r - 0.5*sigma^2)*sqrt(dt))

  # stock price
  S <- matrix(0,ncol=N+1,nrow=N+1) 
  S <- lapply(0:N, function(N.) stock * u^(0:N.) * d^(N.:0))

  # option price
  P <- vector("list", length(S))
  P[[length(P)]] <- ifelse(z*(S[[length(S)]]-strike) > 0, z*(S[[length(S)]]-strike), 0)

  if(american) {
    for(i in N:1)
      P[[i]] <- mapply(max, dis*(p*P[[i+1]][1+(1:i)] + (1-p)*P[[i+1]][1:i]), z*(S[[i]]-strike))
  } else {
    for(i in N:1)
      P[[i]] <- dis*(p*P[[i+1]][1+(1:i)] + (1-p)*P[[i+1]][1:i])
  }
  g <- list()
  g$delta <- diff(P[[2]]) / (S[[1]]*u - S[[1]]*d)
  g$gamma <- (diff(P[[3]][2:3]) / (stock * u^2 - stock * u * d) - diff(P[[3]][1:2]) / (stock*u*d-stock*d^2)) /(0.5*(stock*u^2-stock*d^2))
  g$theta <- (P[[3]][2] - P[[1]]) / (2*dt) / 365
  structure(list(Stock=S,Price=P,Greeks=g), class=c("binomialtree", "JRtree"))
}

LRtree <- function(type=c('ac','ap','ec','ep'),
                   S, #stock=42,
                   X, #strike=40,
                   r, #r=0.1,
                   b, #b=0.1,
                   v, #sigma=0.2,
                   Time, #T=0.5,
                   N=53, 
                   method=2, 
                   force.odd=TRUE, ...) 
{
  stock <- S
  strike <- X
  sigma <- v
  T <- Time

  type <- match.arg(type)
  z <- ifelse(type %in% c('ac','ec'), 1, -1)
  american <- ifelse(type %in% c('ac','ap'), TRUE, FALSE)

  if(force.odd)
    N <- N + (N+1) %% 2  # N should be odd, so we round up if even

  d1 <- ( log(stock / strike) + (b + sigma^2 / 2)*T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)

  if(method==1) {
    # Preizer-Pratt inversion method 1
    hd1 <- 0.5 + sign(d1) * sqrt(0.25-0.25*exp(-(d1/ (N+1/3))^2 * (N+1/6)))
    hd2 <- 0.5 + sign(d2) * sqrt(0.25-0.25*exp(-(d2/ (N+1/3))^2 * (N+1/6)))
  } else if(method==2) {
    # Preizer-Pratt inversion method 2
    hd1 <- 0.5 + sign(d1) * sqrt(0.25-0.25*exp(-(d1/ (N+1/3+0.1/(N+1)))^2 * (N+1/6)))
    hd2 <- 0.5 + sign(d2) * sqrt(0.25-0.25*exp(-(d2/ (N+1/3+0.1/(N+1)))^2 * (N+1/6)))
  } else if(method==3) {
    stop('Method 3 not yet supported')
  #  a <- N-0
  #  b <- 0
  #  (b/a)^2 * ((sqrt((9*a − 1)*(9*b − 1) + 3*z*sqrt(a*(9*b - 1)^2 + b*(9*a - 1)^2 - 9*a*b*z^2))/((9*b-1)^2 - 9*b*z^2))^(1/3))
  }
  dt <- T/N
  p <- hd2
  u <- exp(b*dt) * hd1/hd2
  d <- (exp(b*dt) - p*u) / (1-p)
  dis <- exp(-r*dt)

  # stock price
  S <- matrix(0,ncol=N+1,nrow=N+1) 
  S <- lapply(0:N, function(N.) stock * u^(0:N.) * d^(N.:0))

  # option price
  P <- vector("list", length(S))
  P[[length(P)]] <- ifelse(z*(S[[length(S)]]-strike) > 0, z*(S[[length(S)]]-strike), 0)

  if(american) {
    for(i in N:1)
      P[[i]] <- mapply(max, dis*(p*P[[i+1]][1+(1:i)] + (1-p)*P[[i+1]][1:i]), z*(S[[i]]-strike))
  } else {  # european
    for(i in N:1)
      P[[i]] <- dis*(p*P[[i+1]][1+(1:i)] + (1-p)*P[[i+1]][1:i])
  }
  g <- list()
  g$delta <- diff(P[[2]]) / (stock*u - stock*d)
  g$gamma <- (diff(P[[3]][2:3]) / (stock * u^2 - stock * u * d) - diff(P[[3]][1:2]) / (stock*u*d-stock*d^2)) /(0.5*(stock*u^2-stock*d^2))
  g$theta <- (P[[3]][2] - P[[1]]) / (2*dt) / 365
  structure(list(Stock=S,Price=P,Greeks=g), class=c("binomialtree", "LRtree"))
}


# generic/unused template
BinomialTree <- function(stock=100,strike=100,r=0.06,sigma=0.165,N=3) {
  # I think this one is wrong
  dt <- T/N
  u <- 1.1
  d <- 1/u
  dis <- exp(-r*dt)
  p <- ( exp(r*dt) - d )/(u-d)

  S <- matrix(0,ncol=N+1,nrow=N+1) 
  # stock price
  S <- lapply(0:N, function(N.) stock * u^(0:N.) * d^(N.:0))
  S
  P <- vector("list", length(S))
  P[[length(P)]] <- ifelse(strike-S[[length(S)]] > 0, strike-S[[length(S)]], 0)
  for(i in N:1)
    P[[i]] <- mapply(max, dis*(p*P[[i+1]][1+(1:i)] + (1-p)*P[[i+1]][1:i]), 100-S[[i]])
  list(Stock=S,Put=P)
}

# add in test results from Haug worksheet output as 'correct', use to verify this code
CRRtree_discrete_dividends <- function() {}
JRtree_discrete_dividends <- function() {}
LRtree_discrete_dividends <- function() {}
