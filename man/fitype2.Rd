\name{fitype2}
\alias{fitype2}
\title{Computing the Fisher information matrix under progressive type-II censoring scheme}
\description{Computes the Fisher information matrix under progressive type-I interval censoring scheme. The Fisher information matrix is given by
\deqn{I_{rs}=-E\Bigl(\frac{\partial^2 l(\Theta)}{\partial \theta_r \partial \theta_s}\Bigr),}
where
\deqn{l(\Theta)=\log L(\Theta) \propto C \sum_{i=1}^{m} \log f(x_{i:m:n}{{;}}\Theta) + \sum_{i=1}^{m} R_i \log \bigl[1-F(x_{i:m:n}{{;}}\Theta)\bigr],}
in which \eqn{F(.;\Theta)} is the family cumulative distribution function for \eqn{\Theta=(\theta_1,\dots,\theta_k)^T} and \eqn{r,s=1,\dots,k}, and \eqn{C=n(n-R_1-1)(n-R_1-R_2-2)\dots (n-R_1-R_2-\dots R_{m-1}-m+1)}.
}
\usage{fitype2(plan, param, mle, cdf, pdf, lb = 0, ub = Inf, N = 100)}
\arguments{
\item{plan}{Censoring plan for progressive type-II censoring scheme. It must be given as a \code{data.frame} that includes number of failed items \code{X}, and vector of removed items \code{R}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the maximum likelihood estimators.}
\item{cdf}{Expression for the cumulative distribution function.}
\item{pdf}{Expression for the probability density function.}
\item{lb}{Lower bound of the family support. That is zero by default.}
\item{ub}{Upper bound of the family support. That is \code{Inf} by default.}
\item{N}{An even integer value indicating the number of subdivisions for applying Simpson's integration method.}
}
\details{For some families of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \code{lb} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}
\value{Matrices that represent the expected and observed Fisher information matrices.}
\references{
N. Balakrishnan and AHMED Hossain 2007. Inference for the Type II generalized logistic distribution under progressive Type II censoring, \emph{Journal of Statistical Computation and Simulation}, 77(12), 1013-1031.

M. Teimouri and S. Nadarajah 2016. Bias corrected MLEs under progressive type-II censoring scheme, \emph{Journal of Statistical Computation and Simulation}, 86 (14), 2714-2726.}
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
     N <- 100
  plan <- rtype2(n = n, R = R, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub)
  fitype2(plan = plan, param = param, mle = mle, cdf = cdf, pdf = pdf, lb = lb, ub = ub, N = N)
}
