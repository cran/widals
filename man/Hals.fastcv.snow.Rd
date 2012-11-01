\name{Hals.fastcv.snow}
\alias{Hals.fastcv.snow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ALS Spacial Cross-Validation
}
\description{
Fit Adaptive Least Squares with \eqn{k}-fold cross-validation
}
\usage{
Hals.fastcv.snow(j, rm.ndx, Z, Hs, Ht, Hst.ls, GP.mx)
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
Spacial covariates.  An \eqn{n} x \eqn{ps} numeric matrix.
}
  \item{Ht}{
Temporal covariates.  An \eqn{\tau} x \eqn{pt} numeric matrix.
}
  \item{Hst.ls}{
Space-time covariates.  A list of length \eqn{\tau}, each element containing a \eqn{n} x \eqn{pst} numeric matrix.
}
  \item{GP.mx}{
Hyperparameters.  A \eqn{k.glob} x 2 non-negative matrix.  See \code{\link{MSS.snow}}.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A \eqn{\tau} x \eqn{n} numeric matrix.  The ALS cross-validated predictions of \code{Z}.
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
\code{\link{Hals.snow}}, \code{\link{MSS.snow}}.
}
\examples{

set.seed(99999)

tau <- 70
n.all <- 14

Hs.all <- matrix(rnorm(n.all), nrow=n.all)
Ht <- matrix(rnorm(tau*2), nrow=tau)
Hst.ls.all <- list()
for(i in 1:tau) { Hst.ls.all[[i]] <- matrix(rnorm(n.all*2), nrow=n.all) }

Hst.combined <- list()
for(i in 1:tau) { 
    Hst.combined[[i]] <- cbind( Hs.all, matrix(Ht[i, ], nrow=n.all, 
    ncol=ncol(Ht), byrow=TRUE), Hst.ls.all[[i]] ) 
}

######## use SSsimple to simulate
sssim.obj <- SS.sim.tv( 0.999, Hst.combined, 0.01, diag(1, n.all), tau )



Z.all <- sssim.obj$Z
Z <- Z.all
n <- n.all

Hst.ls <- Hst.ls.all

Hs <- Hs.all

xrho <- 1/10
xreg <- 1/10

GP.mx <- matrix(c(xrho, xreg), nrow=1)

rm.ndx <- create.rm.ndx.ls(n, 10)

Zcv <- Hals.fastcv.snow(j=1, rm.ndx, Z, Hs, Ht, Hst.ls, GP.mx) 



test.rng <- 20:tau

errs.sq <- (Z - Zcv)^2
sqrt( mean(errs.sq[test.rng, ]) )




## The function is currently defined as
function (j, rm.ndx, Z, Hs, Ht, Hst.ls, GP.mx) 
{
    n <- ncol(Z)
    tau <- nrow(Z)
    rho <- GP.mx[j, 1]
    reg <- GP.mx[j, 2]
    Z.hat <- matrix(NA, tau, n)
    for (drop.ndx in rm.ndx) {
        if (!is.null(Hst.ls)) {
            red.Hst.ls <- list()
            Hst0.ls <- list()
            for (i in 1:tau) {
                red.Hst.ls[[i]] <- Hst.ls[[i]][-drop.ndx, , drop = FALSE]
                Hst0.ls[[i]] <- Hst.ls[[i]][drop.ndx, , drop = FALSE]
            }
        }
        else {
            red.Hst.ls <- NULL
            Hst0.ls <- NULL
        }
        if (!is.null(Hs)) {
            red.Hs <- Hs[-drop.ndx, , drop = FALSE]
            Hs0 <- Hs[drop.ndx, , drop = FALSE]
        }
        else {
            red.Hs <- NULL
            Hs0 <- NULL
        }
        Z.hat[, drop.ndx] <- H.als.b(Z[, -drop.ndx, drop = FALSE], 
            Hs = red.Hs, Ht = Ht, Hst.ls = red.Hst.ls, rho, reg, 
            b.lag = 0, Hs0 = Hs0, Ht0 = Ht, Hst0.ls = Hst0.ls)$Z0.hat
    }
    return(Z.hat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
