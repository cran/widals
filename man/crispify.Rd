\name{crispify}
\alias{crispify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Observation-Space Stochastic Correction
}
\description{
Improve observation-space predictions using 'left over' spacial correlation between model residuals
}
\usage{
crispify(locs1, locs2, Z.delta, z.lags.vec, geodesic, alpha, flatten, self.refs, 
lags, stnd.d = FALSE, log10cutoff = -16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locs1}{
Locations of supporting sites.  An \emph{n} x 3 matrix, first column is spacial \eqn{x}, second column is spacial \eqn{y}, third column contains relative temporal 'distance'.  If the \code{geodesic} is \code{TRUE}, make sure latitude is in the first column.
}
  \item{locs2}{
Locations of interpolation sites.  An \eqn{n}* x 3 matrix, where \eqn{n}* is the number of interpolation sites.  See \code{locs1} above.
}
  \item{Z.delta}{
Observed residuals.  A \eqn{\tau} x \eqn{x} matrix.
}
  \item{z.lags.vec}{
Temporal lags.  An integer vector or scalar.
}
  \item{geodesic}{
Use geodesic distance?  Boolean.  If true, distance (used internally) is in units kilometers.
}
  \item{alpha}{
The WIDALS distance rate hyperparameter.  A scalar non-negative number.
}
  \item{flatten}{
The WIDALS 'flattening' hyperparameter.  A scalar non-negative number.  Typically between 0 and some number slightly greater than 1.  When 0, no crispification.
}
  \item{self.refs}{
Which sites are self-referencing?  An integer vector of (zero-based) lag indices, OR a scalar set to \code{-1}.  This argument only has meaning when \code{locs1} is identical to \code{locs2}.  If the \code{lags} argument is, say, 0, then it would be pointless to smooth predictions with existing values.  In this case, we can set \code{self.refs = 0}.  If \code{locs1} is NOT the same as \code{locs2}, then set this argument to \code{-1}.
}
  \item{lags}{
Temporal lags.  An integer vector or scalar.  E.g., if the data's time increment is daily, then \code{lags = c(-1,0,1)} would have \code{crispify} smooth today's predictions using yesterdays, today's, and tomorrow's observed residuals.
}
  \item{stnd.d}{
Spacial compression.  Boolean.
}
  \item{log10cutoff}{
Weight threshold.  A scalar number.  A value of, e.g., -10, will instruct \code{crispify} to ignore weights less than 10^(-10) when smoothing.
}
}
\details{
This function is called inside \code{\link{widals.predict}} and \code{\link{widals.snow}}.  It may be useful for the user in building their own WIDALS model extensions.
}
\value{
A \eqn{\tau} x \eqn{x} matrix.
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
\code{\link{widals.predict}}, \code{\link{widals.snow}}.
}
\examples{

######### here's a itty-bitty example

######### simulate itty-bitty data

tau <- 21 #### number of time points

d.alpha <- 2
R.scale <- 1
sigma2 <- 0.01
F <- 1
Q <- 0

n.all <- 14 ##### number of spacial locations

set.seed(9999)
locs.all <- cbind(runif(n.all, -1, 1), runif(n.all, -1, 1)) #### random location of sensors
D.mx <- distance(locs.all, locs.all, FALSE) #### distance matrix

#### create measurement variance using distance and covariogram
R.all <- exp(-d.alpha*D.mx) + diag(sigma2, n.all)

Hs.all <- matrix(1, n.all, 1) #### constant mean function

##### use SSsimple to simulate system
xsssim <- SS.sim(F=F, H=Hs.all, Q=Q, R=R.all, length.out=tau, beta0=0)
Z.all <- xsssim$Z ###### system observation matrix


######## suppose use the global mean as a prediction

z.mean <- mean(Z.all)

Z.delta <- Z.all - z.mean


z.lags.vec <- rep(0, n.all)

geodesic <- FALSE
alpha <- 5
flatten <- 1

## emmulate cross-validation, i.e., don't use observed site values to predict themselves (zero-based)
self.refs <- 0 
lags <- 0

locs1 <- cbind(locs.all, rep(0, n.all))
locs2 <- cbind(locs.all, rep(0, n.all))

Z.adj <- crispify(locs1, locs2, Z.delta, z.lags.vec, geodesic, alpha, 
    flatten, self.refs, lags, stnd.d = FALSE, log10cutoff = -16) 

Z.adj

Z.hat <- z.mean + Z.adj

sqrt( mean( (Z.all - Z.hat)^2 ) )


######### set flatten to zero -- this means no crispification

Z.adj <- crispify(locs1, locs2, Z.delta, z.lags.vec, geodesic, alpha, 
    flatten=0, self.refs, lags, stnd.d = FALSE, log10cutoff = -16) 

Z.adj

Z.hat <- z.mean + Z.adj

sqrt( mean( (Z.all - Z.hat)^2 ) )






## The function is currently defined as
function (locs1, locs2, Z.delta, z.lags.vec, geodesic, alpha, 
    flatten, self.refs, lags, stnd.d = FALSE, log10cutoff = -16) 
{
    n1 <- nrow(locs1)
    n2 <- nrow(locs2)
    tau <- nrow(Z.delta)
    n.Zd <- ncol(Z.delta)
    Z.out <- rep(0, tau * n1)
    z.rep.in <- rep(0:(n.Zd - 1), length.out = n2)
    t.start <- max(0, -min(lags))
    t.stop <- min(tau, tau - max(lags))
    t.s.s <- c(t.start, t.stop)
    if (geodesic) {
        rlocs1 <- pi * locs1[, 1:2]/180
        rlocs2 <- pi * locs2[, 1:2]/180
        Z.out <- .C("crispify", as.double(cos(rlocs1[, 1]) * 
            cos(rlocs1[, 2])), as.double(cos(rlocs1[, 1]) * sin(rlocs1[, 
            2])), as.double(sin(rlocs1[, 1])), as.double(locs1[, 
            3]), as.double(cos(rlocs2[, 1]) * cos(rlocs2[, 2])), 
            as.double(cos(rlocs2[, 1]) * sin(rlocs2[, 2])), as.double(sin(rlocs2[, 
                1])), as.double(locs2[, 3]), as.double(as.vector(Z.delta)), 
            as.double(Z.out), as.double(alpha), as.double(flatten), 
            as.integer(self.refs), as.integer(length(self.refs)), 
            as.integer(z.lags.vec), as.integer(z.rep.in), as.integer(n.Zd), 
            as.integer(n1), as.integer(n2), as.integer(tau), 
            as.integer(stnd.d), as.integer(t.s.s), as.integer(geodesic), 
            as.double(log10cutoff))[[10]]
    }
    else {
        Z.out <- .C("crispify", as.double(locs1[, 1]), as.double(locs1[, 
            2]), as.double(0), as.double(locs1[, 3]), as.double(locs2[, 
            1]), as.double(locs2[, 2]), as.double(0), as.double(locs2[, 
            3]), as.double(as.vector(Z.delta)), as.double(Z.out), 
            as.double(alpha), as.double(flatten), as.integer(self.refs), 
            as.integer(length(self.refs)), as.integer(z.lags.vec), 
            as.integer(z.rep.in), as.integer(n.Zd), as.integer(n1), 
            as.integer(n2), as.integer(tau), as.integer(stnd.d), 
            as.integer(t.s.s), as.integer(geodesic), as.double(log10cutoff))[[10]]
    }
    dim(Z.out) <- c(tau, n1)
    return(Z.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
