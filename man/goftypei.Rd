\name{goftypei}
\alias{goftypei}
\title{Computing goodness-of-fit (GOF) measures under progressive type-I interval censoring scheme.}
\description{The goodness-of-fit (GOF) measures consist of Anderson-Darling (\code{AD}) and Cram\'{e}er-von Misses (\code{CVM}) statistics for progressive type-I interval censoring scheme are given, respectively, by
\deqn{AD=n\sum_{i=1}^{m}\gamma^{2}_{i}\log\left[\frac{A_{i+1}\bigl(1-A_i\bigr)}{A_i\bigl(1-A_{i+1}\bigr)}\right]+2n\sum_{i=1}^{m}\gamma_{i}\log\Bigl(\frac{1-A_{i+1}}{1-A_i}\Bigr)-n\bigl(A_{m+1}-A_1\bigr)}
\deqn{-n\log\Bigl(\frac{1-A_{m+1}}{1-A_1}\Bigr)+n\bigl(1-A_{m+1}-\log A_{m+1}\bigr),}
\deqn{{CVM}=n\sum_{i=1}^{m}\gamma^{2}_{i}\bigl(A_{i+1}-A_i\bigr)-n\sum_{i=1}^{m}\gamma_{i}\bigl(A^{2}_{i+1}-A^2_i\bigr)+\frac{n}{3}\bigl(A^{3}_{m+1}-A^{3}_{1}\bigr)+\frac{n}{3}\bigl(1-A_{m+1}\bigr)^3,}
where \eqn{R_0=0}, \eqn{\gamma_{i}=\bigl(\sum_{j=1}^{i}{X_j}+\sum_{j=1}^{i-1}{R_j}\bigr)/n}, and \eqn{A_i=G\bigl(T_{i-1}\big|\widehat{\Theta}\bigr)}, for \eqn{i=1,\dots,m}.
}
\usage{goftypei(plan, param, mle, cdf.expression=TRUE, pdf.expression=FALSE, cdf, pdf,
lower = 0)}
\arguments{
\item{plan}{Censoring plan for progressive type-I interval censoring scheme. It must be given as a \code{data.frame} that includes vector of upper bounds of the censoring times \code{T}, vector of number of failed subjects \code{X}, and vector of removed subjects in each interval \code{R}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{Vector of the estimated parameters.}
\item{cdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the cumulative distribution function.}
\item{pdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{lower}{Lower bound of the family support. That is 0, by default.}
}
\details{We note that for lifetime distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \eqn{lower} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf. Theoretically, for lifetime distribution, we have \code{lower}=\eqn{T_{0}=0}.}
\value{A vector of GOF measures consist of Anderson-Darling (\code{AD}) and Cram\'{e}er-von Misses (\code{CVM}) statistics.}
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
pdf <- quote( exp(beta*(1-exp(x^alpha)))*(beta*(exp(x^alpha)*(x^(alpha-1)*alpha))) )
cdf <- quote( 1-exp(beta*(1-exp(x^alpha))) )
param <- c("alpha","beta")
mle <- c(0.4, 0.05)
lower <- 1e-4
goftypei(plan = plan, param = param, mle = mle, cdf.expression = TRUE, pdf.expression = FALSE,
cdf = cdf, pdf = pdf, lower = lower)
}
