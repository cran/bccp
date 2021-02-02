\name{fitype1}
\alias{fitype1}
\title{Computing the Fisher information matrix under progressive type-I interval censoring scheme}
\description{Computes the Fisher information matrix under progressive type-I interval censoring scheme. The Fisher information matrix is given by
\deqn{I_{rs}=-E\Bigl(\frac{\partial^2 l(\Theta)}{\partial \theta_r \partial \theta_s}\Bigr),}
where
\deqn{l(\Theta)=\log L(\Theta) \propto \sum_{i=1}^{m}X_i \log \bigl[F(t_{i}{{;}}\Theta)-F(t_{i-1}{{;}}\Theta)\bigr]+\sum_{i=1}^{m}R_i\bigl[1-F(t_{i}{{;}}\Theta)\bigr],}
in which \eqn{F(.;\Theta)} is the family cumulative distribution function for \eqn{\Theta=(\theta_1,\dots,\theta_k)^T} and \eqn{r,s=1,\dots,k}.
}
\usage{fitype1(plan, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf
        , lb = 0)}
\arguments{
\item{plan}{Censoring plan for progressive type-I interval censoring scheme. It must be given as a \code{data.frame} that includes vector of upper bounds of the censoring times \code{T}, vector of number of failed subjects \code{X}, vector of removed subjects in each interval \code{R}, and percentage of the removed alive items in each interval \code{P}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the maximum likelihood estimators.}
\item{cdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the cumulative distribution function.}
\item{pdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{lb}{Lower bound of the family support. That is zero by default.}
}
\details{For some families of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \code{lb} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}
\value{Matrices that represent the expected and observed Fisher information matrices.}
\references{
N. Balakrishnan and E. Cramer. 2014. \emph{The art of progressive censoring}. New York: Springer.

D. G. Chen and Y. L. Lio 2010. Parameter estimations for generalized exponential distribution under progressive
type-I interval censoring, \emph{Computational Statistics and Data Analysis}, 54, 1581-1591.

M. Teimouri 2020. Bias corrected maximum likelihood estimators under progressive type-I interval censoring scheme, \emph{Communications in Statistics-Simulation and Computation}, doi.org/10.1080/036
10918.2020.1819320

}

\author{Mahdi Teimouri}
\examples{
  data(plasma)
      n <- 20
  param <- c("alpha","beta")
    mle <- c(0.4, 0.05)
    cdf <- quote( 1-exp( beta*(1-exp( x^alpha )) ) )
       pdf <- quote( exp( beta*(1-exp( x^alpha )) )*( beta*(exp( x^alpha )*(x^(alpha-1)*alpha) )) )
     lb <- 0
     plan <- rtype1(n = n, P = plasma$P, T = plasma$upper, param = param, mle = mle, cdf.expression
                  = FALSE, pdf.expression = TRUE, cdf = cdf, pdf = pdf, lb = lb)
fitype1(plan = plan, param = param, mle = mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf =
        cdf, pdf = pdf, lb = lb)
}
