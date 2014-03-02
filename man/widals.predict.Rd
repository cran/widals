\name{widals.predict}
\alias{widals.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WIDALS Interpolation
}
\description{
Interpolate to unmonitored sites using WIDALS
}
\usage{
widals.predict(Z, Hs, Ht, Hst.ls, locs, lags, b.lag, Hs0, Hst0.ls, locs0, 
geodesic = FALSE, wrap.around = NULL, GP, stnd.d = FALSE, ltco = -16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Space-time data.  A \eqn{\tau} x \eqn{n} numeric matrix.
}
  \item{Hs}{
Spacial covariates (of supporting sites).  An \eqn{n} x \eqn{p_s} numeric matrix.
}
  \item{Ht}{
Temporal covariates (of supporting sites).  A \eqn{\tau} x \eqn{p_t} numeric matrix.
}
  \item{Hst.ls}{
Space-time covariates (of supporting sites).  A list of length \eqn{\tau}, each element should be a \eqn{n} x \eqn{p_st} numeric matrix.
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
  \item{Hs0}{
Spacial covariates (of interpolation sites).  An \eqn{n}* x \eqn{p_s} matrix, or \code{NULL}.
}
  \item{Hst0.ls}{
Space-time covariates (of interpolation sites).  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{p_st} matrix.
}
  \item{locs0}{
Locations of interpolation sites.  An \emph{n}* x 2 numeric matrix.  See \code{locs} argument above.
}
  \item{geodesic}{
Use geodesic distance?  Boolean.  If true, distance (used internally) is in units kilometers.
}
  \item{wrap.around}{
**Unused.
}
  \item{GP}{
Widals hyperparameters.  A non-negative vector.
}
  \item{stnd.d}{
Spacial compression.  Boolean.
}
  \item{ltco}{
Weight threshold.  A scalar number.  A value of, e.g., -10, will instruct \code{crispify} to ignore weights less than 10^(-10) when smoothing.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A \eqn{\tau} x \eqn{n}* matrix.  The WIDALS predictions at \code{locs0}.
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
\code{\link{crispify}}, \code{\link{H.als.b}}, \code{\link{widals.snow}}.
}


\examples{
	
#### similar to example provided in H.als.b.
	
set.seed(99999)

library(SSsimple)

tau <- 70
n.all <- 14

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


ndx.support <- 1:10
ndx.interp <- 11:14

locs <- locs.all[ndx.support, ]
locs0 <- locs.all[ndx.interp, ]

Z.all <- sssim.obj$Z
Z <- Z.all[ , ndx.support]
Z0 <- Z.all[ , ndx.interp]

Hst.ls <- subsetsites.Hst.ls(Hst.ls.all, ndx.support)
Hst0.ls <- subsetsites.Hst.ls(Hst.ls.all, ndx.interp)

Hs <- Hs.all[ ndx.support, , drop=FALSE]
Hs0 <- Hs.all[ ndx.interp, , drop=FALSE]

test.rng <- 20:tau


################# use ALS
xrho <- 1/10
xreg <- 1/10
xALS <- H.als.b(Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, rho=xrho, reg=xreg, 
b.lag=-1, Hs0=Hs0, Ht0=Ht, Hst0.ls=Hst0.ls) 

errs.sq <- (Z0 - xALS$Z0.hat)^2
sqrt( mean(errs.sq[test.rng, ]) )

################# now use WIDALS

GP <- c(1/10, 1/10, 2, 0, 1)
Zwid <- widals.predict(Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, locs=locs, lags=c(0), 
b.lag=-1, Hs0=Hs0, Hst0.ls=Hst0.ls, locs0=locs0, FALSE, NULL, GP) 

errs.sq <- (Z0 - Zwid)^2
sqrt( mean(errs.sq[test.rng, ]) )



## The function is currently defined as
function (Z, Hs, Ht, Hst.ls, locs, lags, b.lag, Hs0, Hst0.ls, 
    locs0, geodesic = FALSE, wrap.around = NULL, GP, stnd.d = FALSE, 
    ltco = -16) 
{
    tau <- nrow(Z)
    n <- nrow(locs)
    k <- length(lags)
    n0 <- nrow(locs0)
    rho <- GP[1]
    reg <- GP[2]
    alpha <- GP[3]
    beta <- GP[4]
    flatten <- GP[5]
    locs0.3D <- cbind(locs0, rep(0, n0))
    locs.long.3D <- cbind(rep(locs[, 1], k), rep(locs[, 2], k), 
        beta * rep(lags, each = n))
    z.lags.vec <- rep(lags, each = n)
    ALS <- H.als.b(Z = Z, Hs = Hs, Ht = Ht, Hst.ls = Hst.ls, 
        rho = rho, reg = reg, b.lag = b.lag, Hs0 = Hs0, Ht0 = Ht, 
        Hst0.ls = Hst0.ls)
    Y.als <- ALS$Z.hat
    dim(Y.als)
    Y0.als <- ALS$Z0.hat
    dim(Y0.als)
    rm(ALS)
    Z.delta <- Z - Y.als
    Z.delta <- Z.clean.up(Z.delta)
    Z.adj <- crispify(locs1 = locs0.3D, locs2 = locs.long.3D, 
        Z.delta = Z.delta, z.lags.vec = z.lags.vec, geodesic = geodesic, 
        alpha = alpha, flatten = flatten, self.refs = c(-1), 
        lags = lags, stnd.d = stnd.d, log10cutoff = ltco)
    Z0.wid <- Y0.als + Z.adj
    return(Z0.wid)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
