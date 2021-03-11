black76F <- function(type=c("call","put"),F, X, T, Tf, r, sigma) {
  d1 <- ( log(F/X) + (sigma^2/2)*T ) / (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  if(match.arg(type)=="call")
    exp(-r*Tf) * (F*pnorm(d1) - X*pnorm(d2))
}
