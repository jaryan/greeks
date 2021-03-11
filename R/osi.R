# OSI symbology initiative helper functions

as.osi <- function(x,...) {
  UseMethod("as.osi")
}

as.osi.contract <- function(x, ...) {
  osi <- paste(sprintf("%-6s",x$underlying),format(x$expiry,"%y%m%d"),
               x$right, sprintf("%08s",x$strike*1000), sep="")
  structure(osi, class='osi')
}

as.osi.default <- function(x,...) {
  if(inherits(x,"osi"))
    return(x)
  osi <- paste(sprintf("%-6s",gsub(" ","",substr(sprintf("%21s",x),0,6))),substr(sprintf("%21s",x),7,21),sep="")
  structure(osi, class='osi')
}

is.osi <- function(osi) {
  inherits(osi, "osi")
}

try.osi <- function(osi) {
  try(as.osi(osi), silent=TRUE)
}
# osi print method
print.osi <- function(x, ...) print(unclass(x))

str.osi <- function(object, ...) {
  cat('OSI: ',unclass(object),'\n')
  str(unclass(as.contract(object)))
}

as.data.frame.osi <- function(x, ...) {
  structure(list(x), row.names=.set_row_names(length(x)), class="data.frame")
}


expiry <- function(osi) {
  osi <- try.osi(osi)
  if(!is.osi(osi))
    stop("must be of class 'osi'")
  as.POSIXct(strptime(unclass(substr(osi, 7,12)),"%y%m%d"))
}

right <- function(osi) {
  osi <- try.osi(osi)
  if(!is.osi(osi))
    stop("must be of class 'osi'")
  unclass(substr(osi, 13,13))
}

strike <- function(osi) {
  osi <- try.osi(osi)
  if(!is.osi(osi))
    stop("must be of class 'osi'")
  as.numeric(substr(osi, 14,21)) / 1000
}

tte <- function(osi, from=Sys.Date()) {
  osi <- try.osi(osi)
  if(!is.osi(osi))
    stop("must be of class 'osi'")
  as.integer(as.Date(strptime(unclass(substr(osi, 7,12)),"%y%m%d")) - as.Date(from))
}

underlying <- function(osi) {
  osi <- try.osi(osi)
  if(!is.osi(osi))
    stop("must be of class 'osi'")
  gsub(" ","",unclass(substr(osi, 1,6)))
}

Ops.osi <- function(e1, e2) {
  # construct strategies using contracts/osi

  if(.Generic=='-' && missing(e2)) {
    #e2 <- structure(e2, pos = -1, class = c("osi_list", "osi"))
    e2 <- -1
  }

  if( is.numeric(e1))  {
    structure( e2, pos=e1, class=c("osi_list","osi"))
  } else
  if( is.numeric(e2))  {
    structure( e1, pos=e2, class=c("osi_list","osi"))
  } else
  if(is.list(e1) || is.list(e2)) {
    if( !is.list(e1))
      e1 <- list(e1)
    if( !is.list(e2))
      e2 <- list(e2)
    structure(append(e1,e2), class=c("osi_list","osi"))
  } else {
    structure(list(e1,e2), class=c("osi_list","osi"))
  }
}


