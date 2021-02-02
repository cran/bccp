\name{coxbctype2}
\alias{coxbctype2}
\title{Computing the bias corrected maximum likelihood estimator.}
\description{Computing the bias corrected maximum likelihood estimator (MLE) for the parameters of the general family of distributions under progressive type-II censoring scheme.
}
\usage{coxbctype2(plan, param, mle, cdf, pdf, lb = 0, ub = Inf, N = 100)}
\arguments{
\item{plan}{Censoring plan for progressive type-II censoring scheme. It must be given as a \code{data.frame} that includes number of failed items \code{X}, and vector of removed items \code{R}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the maximum likelihood estimators.}
\item{cdf}{Expression for the cumulative distribution function.}
\item{pdf}{Expression for the probability density function.}
\item{lb}{Lower bound of the family's support. That is zero by default.}
\item{ub}{Upper bound of the family support. That is \code{Inf} by default.}
\item{N}{An even integer value indicating the number of subdivisions for applying Simpson's integration method.}
}
\details{For some families of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \code{lb} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}
\value{A list of the outputs including: a matric that represents the variance-covariance matrix of the MLE, a matrix that represents the variance-covariance matrix of the bias corrected MLE, a list of three outputs including MLE, bias of MLE, and bias corrected MLE, a list of godness-of-fit measures consists of Anderson-Darling (\code{AD}), Cramer-von Misses (\code{CVM}), log-likelihood statistics.}
\references{
D. R. Cox and E. J. Snell 1968. A general definition of residuals. \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 30(2), 248-265.

M. Teimouri and S. Nadarajah 2016. Bias corrected MLEs under progressive type-II censoring scheme, \emph{Journal of Statistical Computation and Simulation}, 86 (14), 2714-2726.
}
\author{Mahdi Teimouri}
\examples{
     n <- 10
     R <- c(5, rep(0, n-6) )
 param <- c("alpha","beta")
   mle <- c(2,6)
   pdf <- quote( alpha/beta*(x/beta)^(alpha-1)*exp( -(x/beta)^alpha ) )
   cdf <- quote( 1-exp( -(x/beta)^alpha ) )
    lb <- 0
    ub <- Inf
     N <- 100
  plan <- rtype2(n = n, R = R, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub)
  coxbctype2(plan = plan, param = param, mle = mle, cdf = cdf, pdf = pdf, lb = lb, ub = ub, N = N)
}
