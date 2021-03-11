\name{tte}
\alias{tte}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Time To Expiration
}
\description{
Find number of calendar days until expiration given a contract specification.
}
\usage{
tte(osi, from = Sys.Date())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{osi}{
A character string formatted per OSI specification or \code{osi} object. 
}
  \item{from}{
A time object or string in CCYY-MM-DD format to be used as reference date. 
}
}
\details{
This calculates the difference between an option contract's expiration date and
the date indicated with \code{from}. 
}
\value{
An integer number of days until expiration. May be negative if \code{from} is before expiration.
}
\references{
OSI Specification
}
\author{
Jeffrey A. Ryan
}
\note{
The value returned is in calendar days, not trading days.
}

\seealso{
\code{\link{osi}}
}
\examples{
tte("AAPL  120122C00450000","2011-12-03")
}
\keyword{ data }
