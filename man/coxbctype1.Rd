 \name{coxbctype1}
\alias{coxbctype1}
\title{Computing the bias corrected maximum likelihood estimator.}
\description{Computing the bias corrected maximum likelihood estimator (MLE) for the parameters of the general family of distributions under progressive type-I interval censoring scheme. Let \eqn{y_1,y_2,\ldots,y_n} represent the lifetimes of \eqn{n} items that independently follow the cumulative distribution function (cdf) \eqn{F(.,\theta_{j})} under a progressive type-I interval censoring scheme. We use \eqn{bctypei} to compute the bias corrected ML estimator using the method of Cox and Snell (1968). Let \eqn{[T_{0}-T_{1})},\eqn{[T_{1}-T_{2})},\ldots,\eqn{[T_{m-1}-T_{m})} show a number of \eqn{m} censoring time intervals, \eqn{{\bf{X}}=(X_{1},X_{2},\ldots,X_{m})} denotes the vector of failed items, and \eqn{{\bf{R}}=(R_{1},R_{2},\ldots,R_{m})} represents the vector of removed alive items in each interval, from \eqn{m\geq 1}. A schematic, given by the following, displays the progressive type-I interval censoring scheme. We note that the sample size \eqn{n} is \eqn{n=\sum_{i=1}^{m}X_{i}+\sum_{i=1}^{m}R_{i}}. Furthermore, \eqn{R_i} can be determined by the pre-specified percentage of the remaining surviving items at \eqn{T_i} or equivalently \eqn{R_i=\lfloor P_i X_i\rfloor}, for \eqn{i=1,\ldots,m}. Here, \eqn{\lfloor z\rfloor} denotes the largest integer less than or equal to \eqn{z}.
\tabular{cccc}{
%\tab number of failures\tab removed alive items\tab\cr
row       \tab Time interval          \tab \eqn{{\bf{X}}}\tab \eqn{{\bf{R}}}\cr
  1       \tab\eqn{[T_{0}, T_{1})}    \tab \eqn{X_1}     \tab      \eqn{R_1}\cr
  2       \tab\eqn{[T_{1}, T_{2})}    \tab \eqn{X_2}     \tab      \eqn{R_2}\cr
         .\tab                       .\tab              .\tab              .\cr
         .\tab                       .\tab              .\tab              .\cr
         .\tab                       .\tab              .\tab              .\cr
\eqn{m-1} \tab\eqn{[T_{m-2}, T_{m-1})}\tab  \eqn{X_{m-1}}\tab  \eqn{R_{m-1}}\cr
\eqn{m}   \tab\eqn{[T_{m-1},T_{m})}   \tab      \eqn{X_m}\tab      \eqn{R_m}\cr
}
}
\usage{coxbctype1(plan, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf
            , lb = 0)}
\arguments{
\item{plan}{Censoring plan for progressive type-I interval censoring scheme. It must be given as a \code{data.frame} that includes vector of upper bounds of the censoring times \code{T}, vector of number of failed items \code{X}, and vector of removed items in each interval \code{R}.}
\item{param}{Vector of the of the family parameter's names.}
\item{mle}{A vector that contains MLE of the parameters.}
\item{cdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the cumulative distribution function.}
\item{pdf.expression}{Logical. That is \code{TRUE}, if there is a closed form expression for the probability density function.}
\item{cdf}{Expression of the cumulative distribution function.}
\item{pdf}{Expression of the probability density function.}
\item{lb}{Lower bound of the family's support. That is zero by default.}
}
\details{For some families of distributions whose support is the positive semi-axis, i.e., \eqn{x>0}, the cumulative distribution function (cdf) may not be differentiable. In this case, the lower bound of the support of random variable, i.e., \code{lb} that is zero by default, must be chosen some positive small value to ensure the differentiability of the cdf.}
\value{A list of the outputs including: a matric that represents the variance-covariance matrix of the MLE, a matrix that represents the variance-covariance matrix of the bias corrected MLE, a list of three outputs including MLE, bias of MLE, and bias corrected MLE, a list of goodness-of-fit measures consists of Anderson-Darling (\code{AD}), Cramer-von Misses (\code{CVM}), and Kolmogorov-Smirnov (\code{KS}) statistics.}
\references{
Z. Chen 2000. A new two-parameter lifetime distribution with bathtub shape or increasing failure rate function, \emph{Statistics & Probability Letters}, 49 (2), 155-161.

D. G. Chen and Y. L. Lio 2010. Parameter estimations for generalized exponential distribution under progressive
type-I interval censoring, \emph{Computational Statistics and Data Analysis}, 54, 1581-1591.

D. R. Cox and E. J. Snell, 1968. A general definition of residuals. \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 30(2), 248-265.

M. Teimouri, 2020. Bias corrected maximum likelihood estimators under progressive type-I interval censoring scheme, https://doi.org/10.1080/03610918.2020.1819320
}
\author{Mahdi Teimouri}
\examples{
 data(plasma, package="bccp")
  plan <- data.frame(T = plasma$upper, X = plasma$X, P = plasma$P, R = plasma$R)
 param <- c("lambda","beta")
   mle <- c(1.4, 0.05)
   pdf <- quote( lambda*(1-exp( -(x*beta)))^(lambda-1)*beta*exp( -(x*beta)) )
   cdf <- quote( (1-exp( -(x*beta)))^lambda )
    lb <- 0
coxbctype1(plan = plan, param = param, mle = mle, cdf.expression = FALSE, pdf.expression = TRUE,
         cdf = cdf, pdf = pdf, lb = lb)
}
