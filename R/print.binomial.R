print.binomialtree <- function(x, tree=FALSE, ...) {
  #> b <- matrix(nc=4,nr=7)
  #> b[4,1] <- CRRtree()$Put[[1]]
  #> b[4+c(-1,1),2] <- CRRtree()$Put[[2]]
  #> b[4+c(-2,0,2),3] <- CRRtree()$Put[[3]]
  #> b[4+c(-3,-1,1,3),4] <- CRRtree()$Put[[4]]
  #  seq(-4+1,by=2,length.out=4)
  op_prices <- x[[2]]
  nc <- length(op_prices)
  p <- matrix(ncol=length(op_prices), nrow=length(op_prices)*2-1)
  p[length(op_prices),1] <- op_prices[[1]]
  for(i in 2:nc)
    p[nc + seq(-i+1, by=2, length.out=i), i] <- op_prices[[i]]
  dimnames(p) <- list(rep('',2*nc-1), rep('',nc))
  print(p, na.print='')
}
