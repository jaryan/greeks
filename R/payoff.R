strike_range <- function(contracts, lower=0.9, upper=1.1) {
  #contracts <- list(...)
  c(min(sapply(contracts, function(X) strike(X) * lower)),
    max(sapply(contracts, function(X) strike(X) * upper))) %/% 1
}
payoff <- function(contracts, premiums, lower=0.9, upper=1.1) {
  xlim <- strike_range(contracts, lower, upper)
  if(!is.list(contracts))
    contracts <- list(contracts)
  cons <- contracts #as.list(...)
  if(missing(premiums))
    premiums <- rep(0,length(cons))
  poff <- lapply(1:length(cons), function(x) {
    premium <- premiums[x]
    x <- as.osi(cons[[x]])
    pos <- attr(x, 'pos')
    if(is.null(pos))
      pos <- 1
    if(pos > 0) { # long
      if(right(x) == "P")
        payoff_longput(x,premium,xlim=xlim) * pos
      else
        payoff_longcall(x,premium,xlim=xlim) * pos
    } else {      # short
      if(right(x) == "P")
        payoff_shortput(x,premium,xlim=xlim) * -pos
      else
        payoff_shortcall(x,premium,xlim=xlim) * -pos
    }
    })
    structure(rowSums(do.call('cbind', poff)), xlim=seq(xlim[1],xlim[2]),
              class="payoff")
}

payoff_longcall <- function(x, premium=0, xlim=strike(x)*c(0.9,1.1)) {
  xlim <- seq(xlim[1], xlim[2])
  if(premium < 0) premium <- abs(premium)
  ifelse(xlim > strike(x), xlim-strike(x)-premium, -premium)
}

payoff_longput <- function(x, premium=0, xlim=strike(x)*c(0.9,1.1)) {
  xlim <- seq(xlim[1], xlim[2])
  if(premium < 0) premium <- abs(premium)
  ifelse(xlim < strike(x), strike(x)-xlim-premium, -premium)
}

payoff_shortcall <- function(x, premium=0, xlim=strike(x)*c(0.9,1.1)) {
  xlim <- seq(xlim[1], xlim[2])
  if(premium < 0) premium <- abs(premium)
  ifelse(xlim > strike(x), strike(x)-xlim+premium, premium)
}

payoff_shortput <- function(x, premium=0, xlim=strike(x)*c(0.9,1.1)) {
  xlim <- seq(xlim[1], xlim[2])
  if(premium < 0) premium <- abs(premium)
  ifelse(xlim < strike(x), xlim-strike(x)+premium, premium)
}
plot.payoff <- function(x, ...) {
  plot(attr(x, 'xlim'), x, 
       ylab='Profit', xlab='Price', type='n', bty='n', axes=F)
  abline(h=0, lty=2, lwd=2, col='red')
  lines(attr(x, 'xlim'), x, col='grey', lwd=3, type='l')
  axis(2, col='grey', las=1)
  axis(1, col='grey', las=1)
}
