\name{nextExp}
\alias{nextExp}
\alias{optionsExp}
\title{
Calculate Next Expiration Date For Derivatives
}
\description{
Used in conjunction with exchange expiration calendars
included in the \pkg{greeks} package to find the
next expiration date given a set of criteria. Currently
only OCC options are available, though this will be user
extensible.
}
\usage{
nextExp(n = 1,
        mindays = 0, 
        from = Sys.Date(), 
        cycle = "monthly", 
        calendar = .OCC,
        expiration = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Which expiration or expirations to return, defaulting to
n=1 meaning the first 'next' expiry from 'from'.
}
  \item{mindays}{
Should there be a minimum number of days until expiration?
Defaults to zero (no minimum).
}
  \item{from}{
Next expiration from which date. Must be of class Date or coercible
to such.
}
  \item{cycle}{
Which expiration cycle should be used. Currently only 'monthly'
is available.
}
  \item{calendar}{
Which product calendar should be used. \sQuote{.OCC} is
the only included calendar at the time of release.
}
  \item{expiration}{
Should the actual expiration be used, e.g. the day following
the last trading day for OCC equity options, or should last
trading date be returned.
}
}
\details{
The initial prototype is designed for US equity option
expirations as given by the OCC calendar. This is a package
variable that is derived from the official OCC sources.
}
\value{
A Date object of length 'n'.
}
\references{
\url{http://www.optionsclearing.com/about/publications/expiration-calendar-2011.jsp}
}
\author{
Jeffrey A. Ryan
}
\note{
It is always possible that dates change for exchanges or individual contracts. Use
your data provider and broker to verify consistency.
}

\examples{
nextExp()

# get next, next cycle
nextExp(n=2)

# change date of calculation
nextExp(from=as.Date("2011-06-02"))

# add minimum number of days to expiry 
nextExp(from=as.Date("2011-06-02"), min=40)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ chron }
\keyword{ ts }
