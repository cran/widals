\name{H.als.b}
\alias{H.als.b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptive Least Squares
}
\description{
Adaptive Least Squares expecially for large spacio-temporal data
}
\usage{
H.als.b(Z, Hs, Ht, Hst.ls, rho, reg, b.lag = -1, Hs0 = NULL, Ht0 = NULL, Hst0.ls = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Space-time data.  A \eqn{\tau} x \eqn{n} numeric matrix.
}
  \item{Hs}{
Spacial covariates (of supporting sites).  An \eqn{n} x \eqn{ps} numeric matrix.
}
  \item{Ht}{
Temporal covariates (of supporting sites).  A \eqn{\tau} x \eqn{pt} numeric matrix.
}
  \item{Hst.ls}{
Space-time covariates (of supporting sites).  A list of length \eqn{\tau}, each element should be a \eqn{n} x \eqn{pst} numeric matrix. 
}
  \item{rho}{
ALS signal-to-noise ratio (SNR).  A non-negative scalar.
}
  \item{reg}{
ALS regularizer.  A non-negative scalar.
}
  \item{b.lag}{
ALS lag.  A scalar integer, typically -1 (\emph{a-prior}), or 0 (\emph{a-posteriori}).
}
  \item{Hs0}{
Spacial covariates (of interpolation sites).  An \eqn{n}* x \eqn{ps} matrix, or \code{NULL}.
}
  \item{Ht0}{
Temporal covariates (of interpolation sites).  A \eqn{\tau} x \eqn{pt} matrix, or \code{NULL}.  If not \code{NULL}, I cannot imagine a scenario where this shouldn't be \code{Ht}.
}
  \item{Hst0.ls}{
Space-time covariates (of interpolation sites).  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{pst} matrix. 
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A named list.
\item{Z.hat}{A \eqn{\tau} x \emph{n} matrix, the \emph{i}th row of which is the ALS prediction of the supporting sites at time \emph{i}.}
\item{B}{A \eqn{\tau} x \eqn{(ps+pt+pst)} matrix, the \emph{i}th row of which is the ALS state (partial slopes) prediction at time \emph{i}.}
\item{Z0.hat}{A \eqn{\tau} x \emph{n}* matrix, the \emph{i}th row of which is the ALS prediction of the interpolation sites at time \emph{i}.}
\item{inv.LHH}{A \eqn{(ps+pt+pst)} x \eqn{(ps+pt+pst)} matrix.  This is the (ALS predicted) covariate precision matrix at time \eqn{\tau}.}
\item{ALS.g}{The ALS gain at time \eqn{\tau}.}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


\examples{

set.seed(99999)

tau <- 70
n.all <- 14

Hs.all <- matrix(rnorm(n.all), nrow=n.all)
Ht <- matrix(rnorm(tau*2), nrow=tau)
Hst.ls.all <- list()
for(i in 1:tau) { Hst.ls.all[[i]] <- matrix(rnorm(n.all*2), nrow=n.all) }

Hst.combined <- list()
for(i in 1:tau) { Hst.combined[[i]] <- cbind( Hs.all, matrix(Ht[i, ], nrow=n.all, ncol=ncol(Ht), byrow=TRUE), Hst.ls.all[[i]] ) }

######## use SSsimple to simulate
sssim.obj <- SS.sim.tv( 0.999, Hst.combined, 0.01, diag(1, n.all), tau )


ndx.support <- 1:10
ndx.interp <- 11:14

Z.all <- sssim.obj$Z
Z <- Z.all[ , ndx.support]
Z0 <- Z.all[ , ndx.interp]

Hst.ls <- subsetsites.Hst.ls(Hst.ls.all, ndx.support)
Hst0.ls <- subsetsites.Hst.ls(Hst.ls.all, ndx.interp)

Hs <- Hs.all[ ndx.support, , drop=FALSE]
Hs0 <- Hs.all[ ndx.interp, , drop=FALSE]

xrho <- 1/10
xreg <- 1/10
xALS <- H.als.b(Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, rho=xrho, reg=xreg, b.lag=-1, Hs0=Hs0, Ht0=Ht, Hst0.ls=Hst0.ls) 



test.rng <- 20:tau

errs.sq <- (Z0 - xALS$Z0.hat)^2
sqrt( mean(errs.sq[test.rng, ]) )


####################### calculate the 'effective standard errors' (actually 'effective prediction errors') of the ALS partial slopes
rmse <- sqrt(mean((Z[test.rng, ] - xALS$Z.hat[test.rng, ])^2))
rmse
als.se <- rmse * sqrt(xALS$ALS.g) * sqrt(diag(xALS$inv.LHH))
cbind(xALS$B[tau, ], als.se, xALS$B[tau, ]/als.se)


## The function is currently defined as
function (Z, Hs, Ht, Hst.ls, rho, reg, b.lag = -1, Hs0 = NULL, 
    Ht0 = NULL, Hst0.ls = NULL) 
{
    tau <- nrow(Z)
    n <- ncol(Z)
    Z.als <- matrix(NA, tau, n)
    if (is.null(Hs)) {
        use.Hs <- FALSE
        d.s <- 0
    }
    else {
        use.Hs <- TRUE
        d.s <- ncol(Hs)
    }
    if (is.null(Ht)) {
        use.Ht <- FALSE
        d.t <- 0
    }
    else {
        use.Ht <- TRUE
        d.t <- ncol(Ht)
    }
    if (is.null(Hst.ls)) {
        use.Hst.ls <- FALSE
        d.st <- 0
    }
    else {
        use.Hst.ls <- TRUE
        d.st <- ncol(Hst.ls[[1]])
    }
    d <- d.s + d.t + d.st
    B <- matrix(0, tau, d)
    reg.mx <- diag(reg, d)
    LHH <- 0
    LHy <- 0
    g <- 1
    use.H0 <- !is.null(Hs0) | !is.null(Ht0) | !is.null(Hst0.ls)
    if (use.H0) {
        if (!is.null(Hs0)) {
            n0 <- nrow(Hs0)
            Z.als.0 <- matrix(NA, tau, n0)
        }
        else {
            n0 <- nrow(Hst0.ls[[1]])
        }
    }
    low.ndx <- max(1, 1 - b.lag)
    top.ndx <- min(tau, tau - b.lag)
    for (i in low.ndx:top.ndx) {
        if (use.Ht) {
            this.Ht.mx <- matrix(Ht[i, ], n, d.t, byrow = TRUE)
        }
        else {
            this.Ht.mx <- NULL
        }
        this.H <- cbind(Hs, this.Ht.mx, Hst.ls[[i]])
        if (use.H0) {
            if (use.Ht) {
                this.Ht0.mx <- matrix(Ht0[i, ], n0, d.t, byrow = TRUE)
            }
            else {
                this.Ht0.mx <- NULL
            }
            this.H0 <- cbind(Hs0, this.Ht0.mx, Hst0.ls[[i]])
        }
        this.HH <- crossprod(this.H)
        this.Hy <- crossprod(this.H, Z[i, ])
        LHH <- LHH + g * (this.HH - LHH)
        LHy <- LHy + g * (this.Hy - LHy)
        inv.LHH <- try(solve(LHH + reg.mx), silent = TRUE)
        if (class(inv.LHH) != "try=error") {
            inv.LHH <- inv.LHH
        }
        else {
            inv.LHH <- matrix(0, d, d)
        }
        B[i, ] <- inv.LHH \%*\% LHy
        Z.als[i, ] <- B[i + b.lag, ] \%*\% t(this.H)
        if (use.H0) {
            Z.als.0[i, ] <- B[i + b.lag, ] \%*\% t(this.H0)
        }
        else {
            Z.als.0 <- NULL
        }
        g <- (g + rho)/(g + rho + 1)
    }
    return(list(Z.hat = Z.als, B = B, Z0.hat = Z.als.0, inv.LHH = inv.LHH, 
        ALS.g = g))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
