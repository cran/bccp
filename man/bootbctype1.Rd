 \name{bootbctype1}
\alias{bootbctype1}
\title{Computing the bias corrected maximum likelihood estimator under progressive type-I interval censoring scheme using the Bootstrap resampling}
\description{Computes the bias corrected maximum likelihood estimator (MLE) under progressive type-I interval censoring scheme using the Bootstrap resampling. It works by obtaining the empirical distribution of the MLE using bootstrap approach and then constructing the percentile confidence intervals (PCI) suggested by DiCiccio and Tibshirani (1987).
}
\usage{bootbctype1(plan, param, mle, cdf, lb = 0, ub = Inf, nboot = 200, coverage = 0.95)}
\arguments{
\item{plan}{Censoring plan for progressive type-I interval censoring scheme. It must be given as a \code{data.frame} that includes vector of upper bounds of the censoring times \code{T}, vector of number of failed subjects \code{X}, vector of removed subjects in each interval \code{R}, and percentage of the removed alive items in each interval \code{P}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the maximum likelihood estimators.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{lb}{Lower bound of the family support. That is zero by default.}
\item{ub}{Upper bound of the family's support. That is \code{Inf} by default.}
\item{nboot}{Number of Bootstrap resampling.}
\item{coverage}{Confidence or coverage level for constructing percentile confidence intervals. That is 0.95 by default.}
}
\value{A list of two parts including: 1- bias of MLE with respect to mean, bias of MLE with respect to median, lower bound of the percentile confidence interval (LPCI), and upper bound of the percentile confidence interval (UPCI) at the given coverage level and 2- covariance matrix of the MLE obtained using bootstraping.}
\references{
T. J. DiCiccio and R. Tibshirani 1987. Bootstrap confidence intervals and bootstrap approximations. \emph{Journal of the American Statistical Association}, 82, 163-170.

A. J. Lemonte, F. Cribari-Neto, and K. L. P. Vasconcellos 2007. Improved statistical inference for the two-parameter Birnbaum-Saunders distribution. \emph{Computational Statistics and Data Analysis}, 51, 4656-4681.
}
\author{Mahdi Teimouri}
\examples{
   data(plasma)
       n <- 112
   param <- c("lambda","beta")
     mle <- c(1.5, 0.05)
     cdf <- quote( (1-exp( -(x*beta)))^lambda )
     pdf <- quote( lambda*(1-exp( -(x*beta)))^(lambda-1)*beta*exp( -(x*beta)) )
      lb <- 0
    plan <- rtype1(n = n, P = plasma$P, T = plasma$upper, param = param, mle = mle,
                  cdf.expression = FALSE, pdf.expression = TRUE, cdf = cdf, pdf = pdf,
                  lb = lb)
      ub <- Inf
   nboot <- 200
coverage <- 0.95
bootbctype1(plan = plan, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub, nboot = nboot,
       coverage = coverage)
}
