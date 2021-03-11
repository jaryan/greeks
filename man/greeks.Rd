\name{greeks}
\alias{greeks}
\title{
Calculate Option Values and Common Risk Sensitivities (Greeks)
}
\description{
Black-Scholes calculation of call and put values, given
standard inputs. Also calculates 14 risk sensitivities,
commonly referred to as \sQuote{Greeks}.
}
\usage{
greeks(S, X, b, r, Time, v)
}
\arguments{
  \item{S}{
Underlying stock price
}
  \item{X}{
Strike price
}
  \item{b}{
Dividend rate
}
  \item{r}{
Risk-free rate
}
  \item{Time}{
Time remaining until expiry
}
  \item{v}{
Estimated volatility
}
}
\details{
Provides a fast and
vectorized way to calculate option prices and all
associated greeks. To accomplish this
the standard Black-Scholes formula is used, along with
formulas for invididual greek values. Both
call and put values are calculated, and returned in a list
of two lists - one for call values and one for put values.

The function returns a list of fifteen values
for each type of contract, calls or puts. These values are
described below.

The estimated \sQuote{value}, via Black-Scholes.

\sQuote{delta} measures the change of an option
value with respect to changes in the underlying price.

\sQuote{gamma} measures the rate of change in delta with respect
to the changes in the underlying price.

\sQuote{vega} measures the sensitivity to volatility changes.

\sQuote{theta} measures the sensitivity to the passage
of time, or \emph{time-decay}.

\sQuote{rho} measures the sensitivity to the applicable
interest rate.

\sQuote{vanna} measures sensitivity of option delta with respect
to changes in volatility.

\sQuote{charm} measures the delta decay - or the change of delta
over time.

\sQuote{zomma} measures gamma change with respect to changes
in volatility.

\sQuote{speed} measures gamma change with respect to changes
in the underlying.

\sQuote{colour} measures gamma decay.

\sQuote{DvegaDtime} measures the change in \sQuote{vega} with respect
to time.

\sQuote{vomma}, also known as vega convexity or volga,
measures vega change with respect to changes in volatility.

\sQuote{dualdelta} is the first derivative of the option price with
respect to strike - also used as the probablility of the option
finishing in the money.
}

\value{
A named list containing two named lists: one for calls and one
for puts.
}
\references{
\url{http://en.wikipedia.org/wiki/Greeks_(finance)}
}
\author{
Jeffrey A. Ryan
}
\note{
The precision of Black-Scholes pricing isn't the reason
for using it, rather it is the speed at which values as well
as approximate greeks can be estimated for a range of inputs.
}


\examples{
greeks(51.03, # underlying price
       55,    # strike
       0,     # dividend rate
       0,     # risk-free rate
       25/360,# time remaining
       0.5)   # volatility estimate
}
\keyword{ models }
