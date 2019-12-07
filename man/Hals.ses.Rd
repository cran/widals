\name{Hals.ses}
\alias{Hals.ses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effective Standard Errors
}
\description{
Calculate the ALS so-called 'effective standard errors'
}
\usage{
Hals.ses(Z, Hs, Ht, Hst.ls, rho, reg, b.lag, test.rng)
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
Space-time covariates (of supporting sites).  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{p_st} matrix. 
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
  \item{test.rng}{
Temporal test range.  A vector of temporal indices of the model test range.
}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A named list.
\item{estimates}{A \eqn{p_s+p_t+p_st} x 2 matrix, each row giving the ALS partial slope estimate/prediction at time \eqn{\tau}, and the 'effective standard error (prediction error)' for the partial slope.}
\item{inv.LHH}{A \eqn{(p_s+p_t+p_st)} x \eqn{(p_s+p_t+p_st)} matrix.  This is the (ALS predicted) covariate precision matrix at time \eqn{\tau}.}
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

## Please see the example in H.als.b

## The function is currently defined as
function (Z, Hs, Ht, Hst.ls, rho, reg, b.lag, test.rng) 
{
    tau <- nrow(Z)
    xALS <- H.als.b(Z = Z, Hs = Hs, Ht = Ht, Hst.ls = Hst.ls, 
        rho = rho, reg = reg, b.lag = b.lag, Hs0 = NULL, Ht0 = NULL, 
        Hst0.ls = NULL)
    rmse <- sqrt(mean((Z[test.rng, ] - xALS$Z.hat[test.rng, ])^2))
    rmse
    als.se <- rmse * sqrt(xALS$ALS.g) * sqrt(diag(xALS$inv.LHH))
    return(list(estimates = cbind(xALS$B[tau, ], als.se), inv.LHH = xALS$inv.LHH, 
        ALS.g = xALS$ALS.g))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
