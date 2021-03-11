\name{LatticePricing}
\alias{CRRtree}
\alias{JRtree}
\alias{LRtree}
\title{
Lattice Based Option Pricing Methods
}
\description{
A collection of lattice-based pricing methods for options.
Currently implementing the Cox-Ross-Rubinstein (1979) [CRRtree],
Jarrow-Rudd (1983) [JRtree], and Leisen-Reimer (1996) [LRtree]
binomial methods.
}
\usage{
CRRtree(type = c("ac", "ap", "ec", "ep"), 
       S, 
       X, 
       r, 
       b, 
       v, 
       Time, 
       N = 52, lambda=1, drift=0, ...)

JRtree(type = c("ac", "ap", "ec", "ep"), 
       S, 
       X, 
       r, 
       b, 
       v, 
       Time, 
       N = 52, ...)

LRtree(type = c("ac", "ap", "ec", "ep"), 
       S, 
       X, 
       r, 
       b, 
       v, 
       Time, 
       N = 53,
       method = 2, force.odd = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
The type of option to be priced. May be one of
\sQuote{ac} (American call), \sQuote{ap} (American put),
\sQuote{ec} (European call), or \sQuote{ep} (European put).
}
  \item{S}{
The stock price.
}
  \item{X}{
The striking price.
}
  \item{r}{
The risk-free rate.
}
  \item{b}{
The dividend rate.
}
  \item{v}{
The volatility of the underlying - e.g. 0.1 is 10%
}
  \item{Time}{
The time to expiration, expressed in fractional years.
}
  \item{N}{
The number of steps in the lattice. For Leisen-Reimer trees
this should be an odd number. Will be enforced if \code{force.odd} is
set as \code{TRUE} for \code{LRtree} (only).
}
  \item{method}{
Which Preizer-Pratt inversion method? 1 or 2.
}
  \item{force.odd}{
Should LR N-steps be forced to odd.
}
  \item{lambda}{
Adjust tilt per Generalized CRR of Chung and Shih (2007). When lambda=0
this is the regular CRR tree.
}
  \item{drift}{
Correct for deep out of the money. When drift=0, this is the standard
CRR tree.
}
  \item{\ldots}{
Unused. 
}
}
\details{
All trees implemented are translations from Espen Haug's excellent
book.

Greeks are calculated from the tree results via finite-difference.
}
\value{

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jeffrey A. Ryan. Based on code from Espen Haug.
}
\note{
Most of these pure \R implementations are as fast as equivelent
code written in C and interfaced to R. The design is meant to be
expended easily for various approaches that one might desire to
implement.
}

\section{Warning}{
As these methods are approximations, the granularity of the computational
lattice determine the accuracy of the method. Literature on the methods
indicate \code{N} must be set sufficiently high to assure convergence. The
tradeoff to higher accuracy comes at a cost of increasing cost of computation.
}

\seealso{
\code{\link{BlackScholesMerton}}
}
\examples{
# Cox-Ross-Rubinstein 
CRRtree(S=42, X=40, r=0.1, b=0.1, v=0.2, Time=0.5, N=5)

# Cox-Ross-Rubinstein with drift
CRRtree(S=42, X=40, r=0.1, b=0.1, v=0.2, Time=0.5, N=5, drift=1)

# Generalized Cox-Ross-Rubinstein
CRRtree(S=42, X=40, r=0.1, b=0.1, v=0.2, Time=0.5, N=5, lambda=2)
CRRtree(S=42, X=40, r=0.1, b=0.1, v=0.2, Time=0.5, N=5, lambda=0.5)

# Jarrow-Rudd
JRtree(S=42, X=40, r=0.1, b=0.1, v=0.2, Time=0.5, N=5)

# Leisen-Reimer
LRtree(S=42, X=40, r=0.1, b=0.1, v=0.2, Time=0.5, N=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
