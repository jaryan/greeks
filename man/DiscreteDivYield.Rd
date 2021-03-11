\name{DiscreteDivYield}
\alias{DiscreteDivYield}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Option Price Given Discrete Dividend Yields Via Recombining Trees
}
\description{
Calculate the price for an American or European option with a series
of discrete dividend yields, using a recombining tree approach as outlined
in Haug (2004).
}
\usage{
DiscreteDivYield(type = c("ac", "ap", "ec", "ep"), 
                 S, 
                 X, 
                 T, 
                 r, 
                 v, 
                 n,
                 d, 
                 dT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{dT}{
%%     ~~Describe \code{dT} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jeffrey A. Ryan adapted from Espen Haug 2004
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

DiscreteDivYield('ac',S=100,X=102,T=0.5,r=0.1,v=0.30,n=10,d=0.15,dT=0.4)

# Example table from 9-8 from Haug 2004 p. 394
\dontrun{
# dividend yields
d <- c(0.01, 0.02, 0.05, 0.10, 0.15, 0.20)
# dividend times
dT <- c(0.1, 0.25, 0.4)

ex9_8 <- matrix(ncol=6,nrow=12)

# calls, sigma=0.15
ex9_8[1:6,1:3] <- sapply(dT, function(dT) sapply(d, function(d) DiscreteDivYield('ac', S=100, X=102, T=0.5, r=0.1, n=500, v=0.15, d=d, dT)$OV[[1]]))
# calls, sigma=0.30
ex9_8[1:6,4:6] <- sapply(dT, function(dT) sapply(d, function(d) DiscreteDivYield('ac', S=100, X=102, T=0.5, r=0.1, n=500, v=0.30, d=d, dT)$OV[[1]]))
# puts, sigma=0.15
ex9_8[7:12,1:3] <- sapply(dT, function(dT) sapply(d, function(d) DiscreteDivYield('ap', S=100, X=102, T=0.5, r=0.1, n=500, v=0.15, d=d, dT)$OV[[1]]))
# puts, sigma=0.30
ex9_8[7:12,4:6] <- sapply(dT, function(dT) sapply(d, function(d) DiscreteDivYield('ap', S=100, X=102, T=0.5, r=0.1, n=500, v=0.30, d=d, dT)$OV[[1]]))
ex9_8
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
