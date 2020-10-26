\name{fitypei}
\alias{fitypei}
\title{Computing the Fisher information matrix under progressive type-I interval censoring scheme}
\description{Computes the Fisher information matrix under progressive type-I interval censoring scheme. The Fisher information matrix is given by
\deqn{I_{rs}=-E\Bigl(\frac{\partial^2 l(\Theta)}{\partial \theta_r \partial \theta_s}\Bigr),}
where
\deqn{l(\Theta)=\log L(\Theta) \propto \sum_{i=1}^{m}X_i \log \bigl[F(t_{i}{{;}}\Theta)-F(t_{i-1}{{;}}\Theta)\bigr]+\sum_{i=1}^{m}R_i\bigl[1-F(t_{i}{{;}}\Theta)\bigr],}
in which \eqn{F(.;\Theta)} is the family cumulative distribution function for \eqn{\Theta=(\theta_1,\dots,\theta_k)^T} and \eqn{r,s=1,\dots,k}.
}
\usage{fitypei(plan, param, mle, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf,
lower = 0)}
\arguments{
\item{plan}{Censoring plan for progressive type-I interval censoring scheme. It must be given as a \code{data.frame} that includes vector of upper bounds of the censoring times \code{T}, vector of number of failed subjects \eqn{X}, vector of removed subjects in each interval \code{R}, and percentage of the removed alive items in each interval \code{P}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the maximum likelihood estimators.}
\item{cdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the cumulative distribution function.}
\item{pdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{lower}{Lower bound of the family support. That is 0, by default.}
}
\details{For some families of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \code{lower} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}
\value{ A matrix that represents the Fisher information matrix.}
\references{
Z. Chen 2000. A new two-parameter lifetime distribution with bathtub shape or increasing failure rate function, \emph{Statistics & Probability Letters}, 49 (2), 155-161.}

\author{Mahdi Teimouri}
\examples{
# bias corrected MLE for Chen distribution with parameters alpha=0.4 and beta=0.05.
T <- c(5.5,10.5,15.5,20.5,25.5,30.5,40.5,50.5,60.5)
X <- c(18, 16, 18, 10, 11, 8, 13, 4,1);
R <- c(1,1,3,0,0,1,2,3,2)
P <- c(0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,1)
plan <- data.frame(T, X, R, P)
# defining the pdf of Chen distribution
pdf <- quote( exp(beta*(1-exp(x^alpha)))*(beta*(exp(x^alpha)*(x^(alpha-1)*alpha))) )
# defining the cdf of Chen distribution
cdf <- quote( 1-exp(beta*(1-exp(x^alpha))) )
param <- c("alpha","beta")
mle <- c(0.4, 0.05)
lower <- 1e-4
fitypei(plan, param, mle, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf,
lower = lower)
}
