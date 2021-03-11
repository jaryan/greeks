# Code adapted from The Complete Guide To Option Pricing Formulas 2nd Edition
# Espen Haug.  Chapter 9. Options On Stocks
#
# Additional code used is for validation/inspiration:
#  http://www.mathfinance.cn/valuation-of-stock-option-with-discrete-dividend/

EscrowedDividend <- function(type=c('call','put'), S, X, b, r, Time, v, d, dT) {
  # Merton 73
  #adj.s <- S-exp(-r*dT)*d
  adj.s <- S - sum(d*exp(-r*dT))
  if(missing(type))
  BlackScholesMerton(,S=adj.s, X=X, r=r, b=b, Time=Time, v=v)
  else
  BlackScholesMerton(type,S=adj.s, X=X, r=r, b=b, Time=Time, v=v)
}
HHD2 <- function (type = c("call", "put"), S, X, b, r, Time, v, d, dT) 
{
    # this isn't really right, but seems ot provide a wickedly good approximation...
    warning("BROKEN: do not use!!!")
    adj.s <- S - sum(exp(-r * dT) * d)
    adj.vol <- v * S/adj.s
    adj.vol <- sqrt(((v^2) * dT + (adj.vol^2) * (Time - dT))/Time)
    if (missing(type)) 
        BlackScholesMerton(, S = adj.s, X = X, r = r, b = b, 
            Time = Time, v = adj.vol)
    else BlackScholesMerton(type, S = adj.s, X = X, r = r, b = b, 
        Time = Time, v = adj.vol)
}

HaugHaugDividend <- function(type=c('c','p'), S, X, b, r, Time, v, d, dT) {
  warning('BROKEN: do not use!!!')
  # Haug, Haug 1998
  adj.s <- S - sum(d*exp(-r*dT))
  adj.vol <- sqrt((S*v/(S - sum(d*exp(-r*dT))))^2 * (dT)+(v^2)*(Time-dT))
  #adj.vol <- sqrt(sum((adj.s*v/(adj.s - sum(d*exp(-r*dT))))^2 * (dT)+(v^2)*(Time-dT)))
  #adj.vol <- adj.vol[length(adj.vol)]
#  adj.vol <- v*S/adj.s
#  adj.vol <- sqrt(((v^2)*dT + (adj.vol^2)*(Time-dT)) / Time)
  if(missing(type))
  BlackScholesMerton(,S=adj.s, X=X, r=r, b=b, Time=Time, v=adj.vol)
  else
  BlackScholesMerton(type,S=adj.s, X=X, r=r, b=b, Time=Time, v=adj.vol)
}

ChrissDividend <- function(type=c('ec','ep'), S, X, b, r, Time, v, d, dT) {
  # Chriss 1997
  adj.s <- S - (d * exp(-r*dT))
  adj.vol <- v*S/adj.s
  if(missing(type))
  BlackScholesMerton(,S=adj.s, X=X, r=r, b=b, Time=Time, v=adj.vol)
  else
  BlackScholesMerton(type=type,S=adj.s, X=X, r=r, b=b, Time=Time, v=adj.vol)
  # ChrissDividend(, 100, 70, 0.06, 0.06, 1, 0.3, 7, 0.5)$call$value
}

BosGairatShepelevaDividend <- function(type=c('ec','ep'), S, X, b, r, Time, v, d, dT) {
  d.sum <- sum(d * exp(-r*dT))
  lns <- log(S)
  lnx <- log( (X+d.sum) * exp(-r * Time) )
  #lnx <- log( (X+d*exp(-r*dT)) * exp(-r*Time))
  z1 <- (lns-lnx) / (v * sqrt(Time))+v*sqrt(Time)/2
  #z2 <- (lns-lnx) / (v * sqrt(Time))+v*sqrt(Time)
  z2 <- z1 + v * sqrt(Time)

  # for i in 1:n
  dTi <- rep(dT, each=length(dT))
  sum1 <- sum(d * exp(-r * dT) * (pnorm(z1) - pnorm(z1-v*dT/sqrt(Time))) )
  # for i in 1:n; for j in 1:n 
  dTj <- rep(dT, length(dT))
  di <- rep(d,each=length(d))
  dj <- rep(d,length(d))
  mint <- pmin(dTi,dTj)
  sum2 <- sum(di*dj * exp(-r*(dTi+dTj)) * (pnorm(z2)-pnorm(z2-2*v*mint)))

  adj.vol <- sqrt(v^2 + v*sqrt(pi/(2*Time)) * (4*exp(z1^2 / 2 - lns) * sum1 + exp(z2^2 / 2 - 2*lns) * sum2))
  adj.s <- S - sum(d*exp(-r*dT))
  as.list(c(unlist(BlackScholesMerton(type,adj.s, X, r, b=b, Time=Time, v=adj.vol)),adj.v=adj.vol))
  # sapply(1:7, function(N) BosGairatShepelevaDividend('c', 100, 100, 0.06, 0.06, N, 0.25, rep(4,N), seq(1,N)-0.5)$call.value)
}

HaugHaugLewisDividend <- function(type=c('ac','ap','ec','ep'),S,X,b,r,Time,v,d,dT) {
  lns <- log(S)
  cp <- ifelse(grepl('c',type), 'c','p')
  american <- ifelse(grepl('a',type), TRUE, FALSE)
  if( !american) {
    exp(-r*dT) * (integrate(function(x) BlackScholesMerton('c',x-d, X,b=b, r, Time-dT,v)$call$value * dlnorm(x,lns+(r-0.5*v^2)*dT,v*sqrt(dT)),d,X+d)$value +
    integrate(function(x) BlackScholesMerton(cp,x-d, X,b=b, r, Time-dT,v)$call$value * dlnorm(x,lns+(r-0.5*v^2)*dT,v*sqrt(dT)),X+d,Inf)$value)
  } else {
    i1 <- function(x) {
            max(x-X, BlackScholesMerton(cp,x-d, X,b=b, r, Time-dT,v)$call$value) * 
                     dlnorm(x,lns+(r-0.5*v^2)*dT,v*sqrt(dT))
          }
    i2 <- function(x) {
            max(x-X, BlackScholesMerton(cp,x-d, X,b=b, r, Time-dT,v)$call$value) * 
                     dlnorm(x,lns+(r-0.5*v^2)*dT,v*sqrt(dT))
          }
    exp(-r*dT) * (integrate(Vectorize(i1), d, X+d)$value+integrate(Vectorize(i2), X+d, 20*S)$value)
  }
  # notes on doing this for multiple dividends:
  #  the process should entail evaluating at t, then t-1, t-2, etc. This is
  #  done by evaluating the integral of the above, n times. HHL suggest
  #  approximating the params of BSM to get the same value of the option
  #  at t (t-1, t-2..) and using that as the function to integrate instead,
  #  thus avoiding n-fold integrals. So far:
  #
  #  1) do the above for the first (last!) dividend
  #  2) adjust the strike - X(adj) = X + d(n) * exp(-r * (dT(n)-Time))
  #     per ex.
  #
  #  > 100 + 4*exp(-0.06*(2.5-3))
  #   [1] 104.1218
  #  > 104.1218 + 4*exp(-0.06*(1.5-3))
  #   [1] 108.4985
  #  > 108.4985 + 4*exp(-0.06*(.5-3))
  #   [1] 113.1458
  #
  #  3) estimate the v(adj) to arrive at the equal price using BSM
  #  4) use this as the input for the second iteration of the integral above
  #     i.e. C(1) = C(BSM) (S, X(adj), v(adj), dT(n-1), Time)
}

DiscreteDividend <- function(type=c('ec','ep','ac','ap'), S, X, b, r, Time, v, d, dT,
                             method=c('BosVandermark','HaugHaugLewis','HaugHaug','Chriss','BosGairatShepeleva','Escrowed'))
{
  method <- match.arg(method)
  do.call(method, list(type=type,S=S,X=X,b=b,r=r,Time=Time,v=v,d=d,dT=dT))
}

BosVandermarkCashDividend <- function(type=c("ec","ep"),S,X,r,b,Time,v,d,dT) {
  # only for european options!
  n <- length(d)
  if(n != length(dT))
    stop(paste(sQuote('d'),"dividends and",sQuote('dT'),"dividend times must",
               "be of equal length")) 

  Xn <- sum( (Time-dT) / Time * d * exp(-r*dT) )
  Xf <- sum( dT / Time * d * exp(-r*dT) )

  if(match.arg(type)=="ec")
    BlackScholesMerton('c',S-Xn,X+Xf * exp(r*Time),r,r,Time,v)$call$value
    #greeks:::call.value(S-Xn, X+Xf * exp(r*Time), r, r, Time, v)$value
  else
    BlackScholesMerton('p',S-Xn,X+Xf*exp(r*Time),r,r,Time,v)$put$value
    #greeks:::put.value(S-Xn, X+Xf * exp(r*Time), r, r, Time, v)$value
}
## ex table 9-6 TCGOPF
## sapply(1:7, function(N) BosVandermarkCashDividend('ec', 100, 100, 0.06, 0.06, N, 0.25, rep(4,N), seq(1,N)-0.5))


BinomialDiscreteDiv <- function(type=c("ac","ap","ec","ep"),S,X,T,r,v,n,d,dT) {
  z <- ifelse(grepl('c',match.arg(type)), 1, -1)
  dt <- T/n
  Df <- exp(-r*dt)
  u <- exp(v*sqrt(dt))
  d <- 1/u
  uu <- u^2
  p <- (exp(r*dt)-d) / (u-d)
  div.amt <- d[1]
  

}

DiscreteDivYield <- function(type=c("ac","ap","ec","ep"),S,X,T,r,v,n,d,dT) {
#  > div  <- c(8, 5, 5, 5, 5) / 100
#  > divT <- c(1, 2, 3, 4, 6) / 10
#  > ss <- DiscreteDivYield("ec",100,101,0.5,.1,.3,100,div,divT)
#  > ss$OV[[1]]
#  [1] 1.051132
#  > ss$St[[1]]
#  [1] 74.93458


  if(length(d) < 1 || is.null(d)) 
    return(CRRtree(type,S,X,T,r,r,v,n))


  dt <- T/n
  Df <- exp(-r*dt)
  U <- exp(v * sqrt(dt))
  D <- 1/U
  UU <- U^2
  p <- (exp(r*dt) - D) / (U-D)

  type <- match.arg(type)
  z <- ifelse(grepl("c",type), 1, -1)
  american <- ifelse(grepl("a",type), TRUE, FALSE)

  i <- 0:(length(d)-1)

  # calculate discrete steps at which dividends must be accounted for
  steps <- as.integer(dT[i+1]/T*n)

  # cumulative dividend adjustments
  sumdiv <- prod(c(1, (1-d[i+1])) )

  St <- vector("list", n+1)
  Ov <- vector("list", n+1)

  # nodes at expiration
  i <- 0:n
  St[[n+1]] <- S*U^i*D^(n-i)*sumdiv
#print(St[[n+1]])
#  St[[n+1]] <- (S*U^i*D^(n-i))-5
#print(St[[n+1]])
  Ov[[n+1]] <- mapply(max, z*(St[[n+1]]-X),0)

  #for(m in length(steps):1) {    ### WORKING on n=6...
  for(m in (n-1):1) {
    if(m %in% steps) {
#      cat('div:',m,'\n')
      div <- d[match(m, steps)]
      St[[m+1]] <- (St[[m+2]]/(1-div)*D)[-1]
#      St[[m+1]] <- ((St[[m+2]]-5)*D)[-1]
      Ov[[m+1]] <- ( p*Ov[[m+2]][-1] + (1-p)*Ov[[m+2]][-length(Ov[[m+2]])]) * Df
    } else {
#      cat('no div:',m,'\n')
      St[[m+1]] <- St[[m+2]][-1]*D
      Ov[[m+1]] <- ( p*Ov[[m+2]][-1] + (1-p)*Ov[[m+2]][-length(Ov[[m+2]])]) * Df
    }
    if(american) {
#print(Ov[[m+1]])
#print(z*(St[[m+1]]-X))
      Ov[[m+1]] <- mapply(max, Ov[[m+1]], z*(St[[m+1]] - X))
    }
  }
  St[[1]] <- St[[2]][-1]*D
  Ov[[1]] <- ( p*Ov[[2]][-1] + (1-p)*Ov[[2]][-length(Ov[[2]])]) * Df
    if(american) {
      Ov[[1]] <- mapply(max, Ov[[1]], z*(St[[1]] - X))
    }
  # TODO greeks
  list(p=p,Df=Df,D=D,steps=steps,St=St,OV=Ov)
}
