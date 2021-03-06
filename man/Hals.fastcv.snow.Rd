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
Spacial covariates.  An \eqn{n} x \eqn{p_s} numeric matrix.
}
  \item{Ht}{
Temporal covariates.  An \eqn{\tau} x \eqn{p_t} numeric matrix.
}
  \item{Hst.ls}{
Space-time covariates.  A list of length \eqn{\tau}, each element containing a \eqn{n} x \eqn{p_st} numeric matrix.
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



library(SSsimple)

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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
