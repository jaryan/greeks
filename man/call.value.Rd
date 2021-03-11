\name{call.value}
\Rdversion{1.1}
\alias{call.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Call Value and Greeks
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
call.value(S, X, b, r, Time, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{Time}{
%%     ~~Describe \code{Time} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(S, X, b, r, Time, v) {
  S <- as.numeric(S)
  X <- as.numeric(X)
  b <- as.numeric(b)
  r <- as.numeric(r)
  Time <- as.numeric(Time)
  v <- as.numeric(v)

  d1 <- (log(S/X) + (r - b + (v^2)/2) * Time)/(v * sqrt(Time))
  d2 <- d1 - v * (sqrt(Time))
  cp <- (S * exp(-b * Time) * pnorm(d1)) - (X * exp(-r * Time) *
        pnorm(d2))
  delta.c <- exp(-b * Time) * pnorm(d1)
  if(identical(cp, NaN)) {
    cp <- delta.c <- 0  
  }
  list(value=cp, delta=delta.c)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
