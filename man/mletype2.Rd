\name{mletype2}
\alias{mletype2}
\title{Computing the maximum likelihood estimator (MLE) for the parameters of the statistical model fitted to a progressive type-II  censoring scheme.}
\description{Computes the MLE of for the parameters of the model fitted to a progressive type-II censoring scheme with likelihood function
\deqn{l(\Theta)=\log L(\Theta) \propto C \sum_{i=1}^{m} \log f(x_{i:m:n}{{;}}\Theta) + \sum_{i=1}^{m} R_i \log \bigl[1-F(x_{i:m:n}{{;}}\Theta)\bigr],}
in which \eqn{F(.;\Theta)} is the family cumulative distribution function for \eqn{\Theta=(\theta_1,\dots,\theta_k)^T} and \eqn{r,s=1,\dots,k}, and \eqn{C=n(n-R_1-1)(n-R_1-R_2-2)\dots (n-R_1-R_2-\dots-R_{m-1}-m+1)}.
}
\usage{mletype2(plan, param, start, cdf, pdf, method = "Nelder-Mead", lb = 0, ub = Inf, N = 100,
        level = 0.05)}
\arguments{
\item{plan}{Censoring plan for progressive type-II censoring scheme. It must be given as a \code{data.frame} including: number of items placed on the test at time zero and a vector that contains number \code{R}, of the removed alive items.}
\item{param}{Vector of the of the family parameter's names.}
\item{start}{Vector of the initial values.}
\item{pdf}{Expression of the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{method}{The method for the numerically optimization that includes one of \code{CG}, \code{Nelder-Mead}, \code{BFGS}, \code{L-BFGS-B}, \code{SANN}.}
\item{lb}{Lower bound of the family's support. That is zero by default.}
\item{ub}{Upper bound of the family's support. That is \code{Inf} by default.}
\item{N}{An even integer number indicating the number of subdivisions for applying Simpson's integration method.}
\item{level}{Significance level for constructing asymptotic confidence interval That is \code{0.05} by default for constructing a \code{95\%} confidence interval.}
}
\value{MLE, standard error of MLE, and asymptotic confidence interval for MLE.}
\references{
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
 level <- 0.05
  plan <- rtype2(n = n, R = R, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub)
mletype2(plan = plan, param = param, start = mle, cdf = cdf, pdf = pdf, method = "Nelder-Mead",
         lb = lb, ub = ub, N = N, level = level)
}
