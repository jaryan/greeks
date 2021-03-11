greeks <-
function(S,X,b,r,Time,v) {
  S <- as.numeric(S)
  X <- as.numeric(X)
  b <- as.numeric(b)
  r <- as.numeric(r)
  Time <- as.numeric(Time)
  v <- as.numeric(v)

  d1 <- (log(S/X) + (r - b + (v^2)/2)*Time ) / (v * sqrt(Time))
  d2 <- d1 - v*(sqrt(Time))
#  cp <- (S * exp(-b*Time) * pnorm(d1)) -
#        (X * exp(-r*Time) * pnorm(d2))
#  pp <- exp(-r*Time) * X * pnorm(-d2) -
#        exp(-b*Time) * S * pnorm(-d1)
  cp <- call.value(S,X,b,r,Time,v)$value
  pp <- put.value(S,X,b,r,Time,v)$value
  delta.c <- exp( (b-r)*Time ) * pnorm(d1)
  delta.p <- -exp( (b-r)*Time ) * pnorm(-d1)
  vega  <- S * exp(-b*Time) * dnorm(d1) * sqrt(Time)
  theta.c <- -exp(-b*Time) * (S*dnorm(d1)*v)/(2*sqrt(Time)) -
           r * X * exp(-r*Time) * pnorm(d2) 
  theta.p <- -exp(-b*Time) * (S*dnorm(d1)*v)/(2*sqrt(Time)) +
           r * X * exp(-r*Time) * pnorm(-d2) 
  rho.c   <- X * Time * exp(-r*Time) * pnorm(d2)
  rho.p   <- -X * Time * exp(-r*Time) * pnorm(-d2)
  gamma <- exp(-b*Time) * ( dnorm(d1)/( S*v*sqrt(Time)) )
  vanna <- -exp(-b*Time) * dnorm(d1) * (d2/v)
  charm.c <- -b*exp(-b*Time)*pnorm(d1) +
             exp(-b*Time)*dnorm(d1)*
               ( (2*(r-b)*Time - d2*v*sqrt(Time))/(2*Time*v*sqrt(Time)) )
  charm.p <- b*exp(-b*Time)*pnorm(-d1) +
             exp(-b*Time)*dnorm(d1)*
               ( (2*(r-b)*Time - d2*v*sqrt(Time))/(2*Time*v*sqrt(Time)) )
  speed <- -(gamma/S) * (d1/(v*sqrt(Time)) + 1)
  zomma <- gamma * ( (d1*d2 - 1)/v )
  colour <- -exp(-b*Time) * (dnorm(d1)/(2*S*Time*v*sqrt(Time))) *
            (2*b*Time + 1 + 
              ((2*(r-b)*Time - d2*v*sqrt(Time))/
               (v*sqrt(Time)))*d1)
  DvegaDtime <- S*exp(-b*Time)*dnorm(d1)*sqrt(Time) *
                (b + ((r-b)*d1)/(v*sqrt(Time)) - ((1+d1*d2)/(2*Time)))
  vomma <- vega * ((d1*d2) / v)
  dualdelta.c  <- -exp(-r*Time) * pnorm(d2)
  dualdelta.p  <- exp(-r*Time) * pnorm(-d2)
  dualgamma <- exp(-r*Time) * ( dnorm(d2)/(X*v*sqrt(Time)) )
  greeks <-
  list(call=list(value=cp,
       delta=delta.c,
       gamma=gamma,
       vega=vega,
       theta=theta.c,
       rho=rho.c,
       vanna=vanna,
       charm=charm.c,
       zomma=zomma,
       speed=speed,
       colour=colour,
       DvegaDtime=DvegaDtime,
       vomma=vomma,
       dualdelta=dualdelta.c,
       dualgamma=dualgamma),

       put=list(value=pp,
       delta=delta.p,
       gamma=gamma,
       vega=vega,
       theta=theta.p,
       rho=rho.p,
       vanna=vanna,
       charm=charm.p,
       zomma=zomma,
       speed=speed,
       colour=colour,
       DvegaDtime=DvegaDtime,
       vomma=vomma,
       dualdelta=dualdelta.p,
       dualgamma=dualgamma)
     )
   class(greeks) <- "greeks"
   greeks
}

BSput <- 
function(S,X,b,r,Time,v) {
  put <- .Call("BSput", S, X, b, r, Time, v, PACKAGE="greeks")
  attr(put, "call") <- match.call()
  put
}

print.greeks <- function(x, ...) {
  str(unclass(x))
}
