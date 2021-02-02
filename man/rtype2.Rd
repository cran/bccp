\name{rtype2}
\alias{rtype2}
\title{Simulating a progressive type-II censoring scheme}
\description{Simulates a progressive type-II censoring scheme.}
\usage{rtype2(n, R, param, mle, cdf, lb = 0, ub = Inf)}
\arguments{
\item{n}{Number of items placed on the test at time zero.}
\item{R}{A vector that contains number \code{R}, of the removed alive items.}
\item{param}{Character vector that contains the name of family's parameters.}
\item{mle}{The maximum likelihood estimators of the parameters.}
\item{cdf}{Expression for the cumulative distribution function.}
\item{lb}{Lower bound of the family's support. That is zero by default.}
\item{ub}{Upper bound of the family's support. That is \code{Inf} by default.}
}
\value{A data frame with two columns including observed liftimes \code{X}, and number of censored items \code{R}.}
\references{
N. Balakrishnan and R. A. Sandhu 1995. A Simple Simulational Algorithm for Generating Progressive Type-II Censored Samples. \emph{The American Statistician}, 49(2), 229-230.
}
\author{Mahdi Teimouri}
\examples{
    n <- 20
    R <- c(9, rep(0, 10) )
param <- c("alpha","beta")
  mle <- c(0.80, 12)
  cdf <- quote( 1-exp( beta*(1-exp( x^alpha )) ) )
   lb <- 0
   ub <- Inf
rtype2(n = n, R = R, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub)
}
