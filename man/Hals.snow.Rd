\name{Hals.snow}
\alias{Hals.snow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit ALS
}
\description{
Fit Adaptive Least Squares
}
\usage{
Hals.snow(j, Z, Hs, Ht, Hst.ls, b.lag, GP.mx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
Index used by \code{\link{snowfall}}.  A scalar integer.  Which row of \code{GP.mx} to use for the ALS hyperparameters, \code{GP}.
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
  \item{b.lag}{
ALS lag.  A scalar integer, typically -1 (\emph{a-prior}), or 0 (\emph{a-posteriori}).
}
  \item{GP.mx}{
Hyperparameters.  A \eqn{k.glob} x 2 non-negative matrix.  See \code{\link{MSS.snow}}.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A \eqn{\tau} x \eqn{n} numeric matrix.  The ALS predictions of \code{Z}.
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
\code{\link{Hals.fastcv.snow}}, \code{\link{MSS.snow}}.
}


\examples{

set.seed(9999)

tau <- 280
n.all <- 35

Hs.all <- matrix(rnorm(n.all), nrow=n.all)
Ht <- matrix(rnorm(tau*2), nrow=tau)
Hst.ls.all <- list()
for(i in 1:tau) { Hst.ls.all[[i]] <- matrix(rnorm(n.all*3), nrow=n.all) }

Hst.combined <- list()
for(i in 1:tau) { Hst.combined[[i]] <- cbind( Hs.all, matrix(Ht[i, ], nrow=n.all, ncol=ncol(Ht), byrow=TRUE), Hst.ls.all[[i]] ) }

######## use SSsimple to simulate
sssim.obj <- SS.sim.tv( 0.999, Hst.combined, 0.1, diag(1, n.all), tau )



Z.all <- sssim.obj$Z
Z <- Z.all
n <- n.all

Hst.ls <- Hst.ls.all

Hs <- Hs.all

xrho <- 1/10
xreg <- 1/10
b.lag <- -1

GP.mx <- matrix(c(xrho, xreg), nrow=1)

Zcv <- Hals.snow(j=1, Z, Hs, Ht, Hst.ls, b.lag, GP.mx) 


test.rng <- 20:tau

errs.sq <- (Z - Zcv)^2
sqrt( mean(errs.sq[test.rng, ]) )


## The function is currently defined as
function (j, Z, Hs, Ht, Hst.ls, b.lag, GP.mx) 
{
    rho <- GP.mx[j, 1]
    reg <- GP.mx[j, 2]
    Z.hat <- H.als.b(Z = Z, Hs = Hs, Ht = Ht, Hst.ls = Hst.ls, 
        rho = rho, reg = reg, b.lag = b.lag, Hs0 = NULL, Ht0 = NULL, 
        Hst0.ls = NULL)$Z.hat
    return(Z.hat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
