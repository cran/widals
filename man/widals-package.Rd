\name{widals-package}
\alias{widals-package}
\alias{widals}
\docType{package}
\title{
Weighting by Inverse Distance with Adaptive Least Squares for Massive Space-Time Data
}
\description{
Fit, forecast, predict massive spacio-temporal data
}
\details{
\tabular{ll}{
Package: \tab widals\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2012-11-08\cr
License: \tab GPL (>=2)\cr
}
The two essential functions are \code{\link{widals.snow}} and \code{\link{widals.predict}}, both contain an Adaptive Least Squares (ALS) prediction stage and complementary 'stochastic adjustment' stage.  The function \code{\link{H.als.b}} solely fits with ALS.

This package offers the user a metaheuristic stochastic search to locate the scalar WIDALS hyperparameters.  The function \code{\link{MSS.snow}} along with helper functions \code{\link{fun.load}} serve this end.  In fairness, providing some useful amount of generality makes this aspect of \code{widals} a bit challenging to learn.  The user new to this package should expect to spend a couple hours playing with the examples before effectively applying these functions to their own data.
}
\author{
Dave Zes

Maintainer: <zesdave@gmail.com>
}

%\references{
%}

\keyword{ package }
\seealso{
%%	\code{\link[<pkg>LatticeKrig<pkg>]{<LatticeKrig>}}
Package \code{LatticeKrig}.
}


%\examples{
%}
