contract <- function(underlying., expiry., strike., right.) {
  # contract("AAPL", "201010", 250, "C")
  if(nchar(underlying.) > 6)
    stop("'underlying' must be 6 characters or less")
  expiry. <- as.Date(expiry.)
  strike. <- as.numeric(strike.)
  right. <- substr(toupper(right.),1,1)
  structure(list(underlying=underlying., expiry=expiry., strike=strike., right=right.),
            class="contract")
}

as.contract <- function(x, ...) {
  UseMethod("as.contract")
}

as.contract.contract <- function(x, ...) {
  x
}

as.contract.osi <- function(x, ...) {
  contract(underlying(x), expiry(x), strike(x), right(x))
}

as.contract.character <- function(x, ...) {
  x <- try.osi(x)
  if(is.osi(x))
    as.contract(x)
  else stop("improperly formatted 'osi' string")
}

as.combo <- function(...) {
  contracts <- lapply(list(...), as.contract)
  structure(contracts, class="combo")
}

is.combo <- function(x) {
  inherits(x, "combo")
}
