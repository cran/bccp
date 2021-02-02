\name{mletype1}
\alias{mletype1}
\title{Computing the maximum likelihood estimator (MLE) for the parameters of the statistical model fitted to a progressive type-I interval censoring scheme.}
\description{Computes the MLE of for the parameters of the model fitted to a progressive type-I interval censoring scheme with likelihood function
\deqn{l(\Theta)=\log L(\Theta) \propto \sum_{i=1}^{m}X_i \log \bigl[F(t_{i}{{;}}\Theta)-F(t_{i-1}{{;}}\Theta)\bigr]+\sum_{i=1}^{m}R_i\bigl[1-F(t_{i}{{;}}\Theta)\bigr],}
in which \eqn{F(.;\Theta)} is the family cumulative distribution function for \eqn{\Theta=(\theta_1,\dots,\theta_k)^T} provided that
\eqn{F(t_{0};\Theta)=0}.
}
\usage{mletype1(plan, param, start, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf,
         method = "Nelder-Mead", lb = 0, ub = Inf, level = 0.05)}
\arguments{
\item{plan}{Censoring plan for progressive type-I interval censoring scheme. It must be given as a \code{data.frame} that includes vector of upper bounds of the censoring times \code{T}, vector of number of failed items \code{X}, and vector of removed items in each interval \code{R}.}
\item{param}{Vector of the of the family parameter's names.}
\item{start}{Vector of the initial values.}
\item{cdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the cumulative distribution function.}
\item{pdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{method}{The method for the numerically optimization that includes one of \code{CG}, \code{Nelder-Mead}, \code{BFGS}, \code{L-BFGS-B}, \code{SANN}.}
\item{lb}{Lower bound of the family's support. That is zero by default.}
\item{ub}{Upper bound of the family's support. That is \code{Inf} by default.}
\item{level}{Significance level for constructing asymptotic confidence interval That is \code{0.05} by default for constructing a \code{95\%} confidence interval.}
}
\value{MLE, standard error of MLE, and asymptotic confidence interval for MLE.}
%\references{}
\author{Mahdi Teimouri}
\examples{
 data(plasma, package="bccp")
   plan <- data.frame(T = plasma$upper, X = plasma$X, P = plasma$P, R = plasma$R)
  param <- c("lambda","beta")
    mle <- c(1.4, 0.05)
    pdf <- quote( lambda*(1-exp( -(x*beta)))^(lambda-1)*beta*exp( -(x*beta)) )
    cdf <- quote( (1-exp( -(x*beta)))^lambda )
     lb <- 0
     ub <- Inf
  level <- 0.05
mletype1(plan = plan, param = param, start = mle, cdf.expression = FALSE, pdf.expression = TRUE,
         cdf = cdf, pdf = pdf, method = "Nelder-Mead", lb = lb, ub = ub, level = level)
}
