\name{widals.snow}
\alias{widals.snow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit WIDALS
}
\description{
Locate the WIDALS hyperparameters
}
\usage{
widals.snow(j, rm.ndx, Z, Hs, Ht, Hst.ls, locs, lags, b.lag, cv = 0, 
geodesic = FALSE, wrap.around = NULL, GP.mx, stnd.d = FALSE, ltco = -16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
Index used by \code{\link{snowfall}}.  A scalar integer.  Which row of \code{GP.mx} to use for the ALS hyperparameters, \code{GP}.
}
  \item{rm.ndx}{
A list of vectors of indices to remove for \emph{k}-fold cross-validation.
}
  \item{Z}{
Data.  A \eqn{\tau} x \eqn{n} numeric matrix.
}
  \item{Hs}{
Spacial covariates.  An \eqn{n} x \eqn{p_s} numeric matrix.
}
  \item{Ht}{
Temporal covariates.  A \eqn{\tau} x \eqn{p_t} numeric matrix.
}
  \item{Hst.ls}{
Space-time covariates.  A list of length \eqn{\tau}, each element containing a \eqn{n} x \eqn{p_st} numeric matrix.
}
  \item{locs}{
Locations of supporting sites.  An \emph{n} x 2 numeric matrix, first column is spacial \eqn{x}, second column is spacial \eqn{y}.  If the \code{geodesic} is \code{TRUE}, make sure latitude is in the first column.
}
  \item{lags}{
Temporal lags for stochastic smoothing.  An integer vector or scalar.  E.g., if the data's time increment is daily, then \code{lags = c(-1,0,1)} would tell the enclosed function \code{\link{crispify}} smooth today's predictions using yesterdays, today's, and tomorrow's observed residuals.
}
  \item{b.lag}{
ALS lag.  A scalar integer, typically -1 (\emph{a-prior}), or 0 (\emph{a-posteriori}).
}
  \item{cv}{
Cross-validation switch.  Currently takes on a value of \code{-2} or \code{2}.  See Details below.
}
  \item{geodesic}{
Use geodesic distance?  Boolean.  If true, distance (used internally) is in units kilometers.
}
  \item{wrap.around}{
**Unused.
}
  \item{GP.mx}{
Hyperparameters.  A \eqn{k.glob} x 2 non-negative matrix.  See \code{\link{MSS.snow}}.
}
  \item{stnd.d}{
Spacial compression.  Boolean.
}
  \item{ltco}{
Weight threshold.  A scalar number.  A value of, e.g., -10, will instruct \code{crispify} to ignore weights less than 10^(-10) when smoothing.
}
}
\details{
When the \code{cv} is set to 2, then this function uses spacial \emph{k}-fold validation, according to the site indices present in \code{rm.ndx}.  When \code{cv} is set to -2, self-referencing sites are given zero-weight, i.e., a site's value is not allowed to contribute to its predicted value.
}
\value{
A \eqn{\tau} x \eqn{n} matrix.  The WIDALS predictions at \code{locs}.
}


%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{crispify}}, \code{\link{H.als.b}}, \code{\link{widals.predict}}.
}
\examples{
	
set.seed(99999)

library(SSsimple)

tau <- 100
n.all <- 35

Hs.all <- matrix(rnorm(n.all), nrow=n.all)
Ht <- matrix(rnorm(tau*2), nrow=tau)
Hst.ls.all <- list()
for(i in 1:tau) { Hst.ls.all[[i]] <- matrix(rnorm(n.all*2), nrow=n.all) }

Hst.combined <- list()
for(i in 1:tau) { 
    Hst.combined[[i]] <- cbind( Hs.all, matrix(Ht[i, ], nrow=n.all, ncol=ncol(Ht), 
    byrow=TRUE), Hst.ls.all[[i]] ) 
}

locs.all <- cbind(runif(n.all, -1, 1), runif(n.all, -1, 1))
D.mx.all <- distance(locs.all, locs.all, FALSE)
R.all <- exp(-2*D.mx.all) + diag(0.01, n.all)

######## use SSsimple to simulate
sssim.obj <- SS.sim.tv( 0.999, Hst.combined, 0.01, R.all, tau )


n <- n.all
locs <- locs.all

Z.all <- sssim.obj$Z
Z <- Z.all


Hst.ls <- Hst.ls.all
Hs <- Hs.all

test.rng <- 20:tau

################  WIDALS, true cross-validation

rm.ndx <- create.rm.ndx.ls(n, 10)

cv <- 2
lags <- c(0)
b.lag <- 0

GP <- c(1/8, 1/12, 5, 0, 1)
GP.mx <- matrix(GP, ncol=length(GP))
Zwid <- widals.snow(j=1, rm.ndx, Z, Hs, Ht, Hst.ls, locs, lags, b.lag, cv = cv, 
geodesic = FALSE, wrap.around = NULL, GP.mx, stnd.d = FALSE, ltco = -16) 

errs.sq <- (Z - Zwid)^2
sqrt( mean(errs.sq[test.rng, ]) )


################  WIDALS, pseudo cross-validation

rm.ndx <- I(1:n)

cv <- -2
lags <- c(0)
b.lag <- -1

GP <- c(1/8, 1/12, 5, 0, 1)
GP.mx <- matrix(GP, ncol=length(GP))
Zwid <- widals.snow(j=1, rm.ndx, Z, Hs, Ht, Hst.ls, locs, lags, b.lag, cv = cv, 
geodesic = FALSE, wrap.around = NULL, GP.mx, stnd.d = FALSE, ltco = -16) 

errs.sq <- (Z - Zwid)^2
sqrt( mean(errs.sq[test.rng, ]) )


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
