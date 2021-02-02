\name{simpson}
\alias{simpson}
\title{Computing integration numerically through the Simpson's method}
\description{Computes the integration for a real-valued function.}
\usage{simpson(fun, lb, ub, N = 100)}
\arguments{
\item{fun}{Integrand expression}
\item{lb}{Lower bound of integration.}
\item{ub}{Upper bound of integration.}
\item{N}{An even integer value indicating the number of subdivisions for applying Simpson's integration method.}
}
\references{
E. Suli and D. Mayers 2003. \emph{An Introduction to Numerical Analysis}, Cambridge University Press.
}
\author{Mahdi Teimouri}
\examples{
    fun <- function(x) 1/sqrt( 2*pi*sigma^2 )*exp( -.5*(x-mu)^2/sigma^2 )
    mu <- 0
 sigma <- 1
    lb <- 0
    ub <- Inf
     N <- 100
simpson(fun = fun, lb = lb, ub = ub, N = N)
}
