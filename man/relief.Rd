\name{relief}
\docType{data}
\alias{relief}
\title{Wingo's pain relief data}
\description{
In order to investigate the effectiveness of an anesthetic
antibiotic ointment in relieving pain caused by wounds, a number of 30 patients take part in a
test where 10 of them are removed from the test under censoring plan \eqn{R=(5,1,4,0)} at times
\eqn{T_{1}=0.25, T_{2}=0.50, T_{3}=0.75} and \eqn{T_{4}=3.50}. We note that the termination time is not fixed, but
it is assumed that \eqn{T_{4}=3.50} sufficiently large to ensure that no observation is right censored,
see Balakrishnan and Cramer (2014).
}
\usage{data(relief)}
\format{A text file with 20 observations.}
\encoding{latin1}
\references{
N. Balakrishnan and E. Cramer 2014. \emph{The Art of Progressive Censoring}. New York, Springer.
}
\keyword{datasets}
\examples{
data(relief)
}
