\name{contract}
\alias{contract}
\alias{as.contract}
\alias{as.contract.osi}
\alias{as.contract.character}
\title{
Convert To And From Class \code{contract}
}
\description{
Constructor function to create an option contract object from
component elements, or from related classes.
}
\usage{
contract(underlying., expiry., strike., right.)

as.contract(x, ...)

\method{as.contract}{character}(x, ...)
\method{as.contract}{osi}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{underlying.}{
The character string of the underlying symbol.
}
  \item{expiry.}{
The expiration date of the contract.  May be a time object
of class POSIXt or Date, or a character string coercible to.
}
  \item{strike.}{
Numeric value up to the limit of the OSI.
}
  \item{right.}{
The right confered, \sQuote{p} or \sQuote{put} for puts,
or \sQuote{c} or \sQuote{call} for calls. Case insensitive.
}
  \item{x}{
  Object to convert.
}
  \item{\ldots}{
  Unused.
} 
}
\details{
This function allows for automatic generation of arbitrary
contract details as required by the OSI identifier format.
There are no checks for sanity beyond the input type validation.
}
\value{
An object of type \code{contract}.
}
\references{
\url{http://www.theocc.com/components/docs/initiatives/symbology/symbology_initiative_v1_8.pdf}
}
\author{
Jeffrey A. Ryan
}


\seealso{
\code{\link{osi}}
}
\examples{
ppl <- contract("PPL",        # Pennsylvania Power and Light underlying
                "2012-10-20", # Expiration date (third Saturday)
                29,           # Striking price
                "c")          # Call
ppl

as.contract("PPL  121020P00029000")

as.contract(as.osi("PPL  121020P00029000"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ data }
