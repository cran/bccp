\name{rtypei}
\alias{rtypei}
\title{Simulating a progressive type-I interval censoring scheme}
\description{Simulates progressive type-I interval censoring scheme when censoring times, number of failed subjects in each interval, number of removed subjects in end of each interval, and percentage of removed items in end of each interval are given in a \eqn{data.frame} structure.}
\usage{rtypei(n, plan, param, mle, cdf.expression = TRUE, pdf.expression = FALSE, cdf,
pdf, lower = 0)}
\arguments{
\item{n}{Number of simulated data.}
\item{plan}{A data frame structure for specifying a progressive type-I interval censoring scheme. It contains vector of upper bounds of the censoring              times \eqn{T}, and percentage of the removed alive items in each interval \eqn{P}.}
\item{param}{A character vector that contains the name of family's parameters.}
\item{mle}{A vector that contains of the maximum likelihood estimators of the parameters.}
\item{cdf.expression}{Logical. That is \eqn{TRUE}, if there is a closed form expression for the cumulative distribution function.}
\item{pdf.expression}{Logical. That is \eqn{TRUE}, if there is a closed form expression for the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{lower}{Lower bound of the family support. That is 0, by default.}
}
\details{For some families of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \eqn{lower} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}
\value{
A data frame with four columns including censoring times, number of failed subjects in each interval, number of removed subjects in end of each interval, and percentage of removed items in end of each interval.
}
\author{Mahdi Teimouri}
\examples{
# bias corrected MLE for Chen distribution with parameters alpha=0.4 and beta=0.05
T <- c(5.5,10.5,15.5,20.5,25.5,30.5,40.5,50.5,60.5)
P <- c(0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,1)
plan <- data.frame(T, P)
pdf <- quote( exp(beta*(1-exp(x^alpha)))*(beta*(exp(x^alpha)*(x^(alpha-1)*alpha))) )
cdf <- quote( 1-exp(beta*(1-exp(x^alpha))) )
param <- c("alpha","beta")
mle <- c(0.4, 0.05)
lower <- 1e-4
n <- 100
rtypei(n, plan, param, mle, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf,
lower = lower)
}
