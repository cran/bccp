\name{bootbctype2}
\alias{bootbctype2}
\title{Computing the bias corrected maximum likelihood estimator under progressive type-I interval censoring scheme using the Bootstrap resampling}
\description{Computes the bias corrected maximum likelihood estimator under progressive type-I interval censoring scheme using the Bootstrap resampling. It works by obtaining the empirical distribution of the MLE using bootstrap approach and then constructing the percentile confidence intervals (PCI) suggested by DiCiccio and Tibshirani (1987).
}
\usage{bootbctype2(plan, param, mle, cdf, pdf, lb = 0, ub = Inf, nboot = 200, coverage = 0.95)}
\arguments{
\item{plan}{Censoring plan for progressive type-II censoring scheme. It must be given as a \code{data.frame} that includes number of failed items \code{X}, and vector of removed items \code{R}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the maximum likelihood estimators.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression for the probability density function.}
\item{lb}{Lower bound of the family support. That is zero by default.}
\item{ub}{Upper bound of the family's support. That is \code{Inf} by default.}
\item{nboot}{Number of Bootstrap resampling.}
\item{coverage}{Confidence or coverage level for constructing percentile confidence intervals. That is 0.95 by default.}
}
\details{For some families of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \code{lb} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}

\value{A list of the outputs including a matrix that represents the variance-covariance matrix of the uncorrected MLE, a matrix that represents the variance-covariance matrix of the corrected MLE, the lower \code{LPCI}, and upped \code{UPCI}, bounds of \code{95\%} percentile confidence interval for \code{param}, the ML estimator, bias value, and bias-corrected estimator. Finally, the goodness-of-fit measures consists of Anderson-Darling (\code{AD}), Cramer-von Misses (\code{CVM}), and Kolmogorov-Smirnov (\code{KS}) statistics.}

\references{
T. J. DiCiccio and R. Tibshirani 1987. Bootstrap confidence intervals and bootstrap approximations. \emph{Journal of the American Statistical Association}, 82, 163-170.

A. J. Lemonte, F. Cribari-Neto, and K. L. P. Vasconcellos 2007. Improved statistical inference for the two-parameter Birnbaum-Saunders distribution.
 \emph{Computational Statistics and Data Analysis}, 51, 4656-4681.
}
\author{Mahdi Teimouri}
\examples{
       n <- 20
       R <- c(9, rep(0, 10) )
   param <- c("alpha","beta")
     mle <- c(0.80, 12)
     cdf <- quote( 1-exp( -(x/beta)^alpha ) )
     pdf <- quote( alpha/beta*(x/beta)^(alpha-1)*exp( -(x/beta)^alpha ) )
      lb <- 0
      ub <- Inf
   nboot <- 200
coverage <- 0.95
    plan <- rtype2(n = n, R = R, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub)
  bootbctype2(plan = plan, param = param, mle = mle, cdf = cdf, pdf = pdf, lb = lb, ub = ub,
          nboot = nboot, coverage = coverage)
}
