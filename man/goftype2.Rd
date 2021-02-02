\name{goftype2}
\alias{goftype2}
\title{Computing goodness-of-fit (GOF) measures under progressive type-II censoring scheme.}
\description{The goodness-of-fit (GOF) measures consist of Anderson-Darling (\code{AD}), Cramer-von Misses (\code{CVM}), and log-likelihood statistics for progressive type-II censoring scheme.}
%are given, respectively, by
%\deqn{AD=n\sum_{i=1}^{m}\gamma^{2}_{i}\log\left[\frac{A_{i+1}\bigl(1-A_i\bigr)}{A_i\bigl(1-A_{i+1}\bigr)}\right]+2n\sum_{i=1}^{m}\gamma_{i}\log\Bigl(\frac{1-A_{i+1}}{1-A_i}\Bigr)-n\bigl(A_{m+1}-A_1\bigr)}
%\deqn{-n\log\Bigl(\frac{1-A_{m+1}}{1-A_1}\Bigr)+n\bigl(1-A_{m+1}-\log A_{m+1}\bigr),}
%\deqn{{CVM}=n\sum_{i=1}^{m}\gamma^{2}_{i}\bigl(A_{i+1}-A_i\bigr)-n\sum_{i=1}^{m}\gamma_{i}\bigl(A^{2}_{i+1}-A^2_i\bigr)+\frac{n}{3}\bigl(A^{3}_{m+1}-A^{3}_{1}\bigr)+\frac{n}{3}\bigl(1-A_{m+1}\bigr)^3,}
%where \eqn{R_0=0}, \eqn{\gamma_{i}=\bigl(\sum_{j=1}^{i}{X_j}+\sum_{j=1}^{i-1}{R_j}\bigr)/n}, and \eqn{A_i=G\bigl(T_{i-1}\big|\widehat{\Theta}\bigr)}, for \eqn{i=1,\dots,m}.
%}
\usage{goftype2(plan, param, mle, cdf, pdf)}

\arguments{
\item{plan}{Censoring plan for progressive type-II censoring scheme. It must be given as a \code{data.frame} that includes number of failed items \code{X}, and vector of removed items \code{R}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the maximum likelihood estimators.}
\item{cdf}{Expression for the cumulative distribution function.}
\item{pdf}{Expression for the probability density function.}
}
\references{
R. Pakyari and N. Balakrishnan 2012. A general purpose approximate goodness-of-fit for progressively Type-II censored
data, \emph{IEEE Transaction on Reliability}, 61, 238-244.
}
\author{Mahdi Teimouri}
\examples{
     n <- 20
     R <- c(5, rep(0, n-6) )
 param <- c("alpha","beta")
   mle <- c(2,6)
   pdf <- quote( alpha/beta*(x/beta)^(alpha-1)*exp( -(x/beta)^alpha ) )
   cdf <- quote( 1-exp( -(x/beta)^alpha ) )
    lb <- 0
    ub <- Inf
  plan <- rtype2(n = n, R = R, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub)
goftype2(plan = plan, param = param, mle = mle, cdf = cdf, pdf = pdf)
}
