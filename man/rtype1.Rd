\name{rtype1}
\alias{rtype1}
\title{Simulating a progressive type-I interval censoring scheme}
\description{Simulates a progressive type-I interval censoring scheme when censoring times and percentage of removed items in end of each interval are given in a \code{data.frame} structure.}
\usage{rtype1(n, P, T, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf,
          lb = 0)}
\arguments{
\item{n}{Number of items placed on the test at time zero.}
\item{P}{A vector of upper bounds of the censoring times \code{T}.}
\item{T}{A vector that contains percentage \code{P}, of the removed alive items in each interval.}
\item{param}{A character vector that contains the name of family's parameters.}
\item{mle}{A vector that contains of the maximum likelihood estimators of the parameters.}
\item{cdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the cumulative distribution function.}
\item{pdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{lb}{Lower bound of the family's support. That is zero by default.}
}
\details{For a family of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \code{lb} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}
\value{
A data frame with four columns including censoring times, number of failed items in each interval, number of removed subjects in end of each interval, and percentage of removed items in end of each interval.
}

\references{
R. Aggarwala 2001. Progressive interval censoring: some mathematical results with applications to inference. \emph{Communications in Statistics-Theory and Methods}, 30(8&9), 1921--1935.
}
\author{Mahdi Teimouri}
\examples{
   data(plasma)
      n <- 112
  param <- c("alpha","beta")
    mle <- c(0.4, 0.05)
    cdf <- quote( 1-exp( beta*(1-exp( x^alpha )) ) )
      pdf <- quote( exp( beta*(1-exp( x^alpha )) )*( beta*(exp( x^alpha )*( x^(alpha-1)*alpha ) )) )
     lb <- 0
rtype1(n = n, P = plasma$P, T = plasma$upper, param = param, mle = mle, cdf.expression = FALSE,
      pdf.expression = TRUE, cdf = cdf, pdf = pdf, lb = lb)
}
