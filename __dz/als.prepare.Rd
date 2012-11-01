\name{als.prepare}
\alias{als.prepare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align Data
}
\description{
Flatten, or reflatten data for lagged auto-regressive ALS
}
\usage{
als.prepare(X, Z, lags, Z.na, tt.rng, rgr = NULL, assignToEnvironment = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Covariates.
}
  \item{Z}{
Responses.
}
  \item{lags}{
Temporal lags.
}
  \item{Z.na}{
Indicator of missing values in \code{Z}.
}
  \item{tt.rng}{
Any temporal index.
}
  \item{rgr}{
Extra inputs.
}
  \item{assignToEnvironment}{
Environment to assign output.  If \code{NULL}, function returns objects as list.
}
}
\details{
This function is somewhat vestigial.  It has no meaningful role in WIDALS, but can be helpful aligning data for the sake of auto-regressive ALS.
}
\value{
If \code{assignToEnvironment} is a valid Environment, then nothing.
If \code{assignToEnvironment} is \code{NULL}, then a named list.
\item{X.sub}{Covariates aligned as prescribed by \code{lags}.}
\item{Y.sub}{Responses aligned as prescribed by \code{lags}.}
\item{Y.na}{\code{Z.na} aligned as prescribed by \code{lags}.}
\item{y.tt.rng}{\code{tt.rng} aligned as prescribed by \code{lags}.}
\item{y.sub.ndx}{The indices that map the rows of \code{Z} to the rows of \code{Y.sub}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ar.align}}, \code{\link{xy.align}}, \code{\link{y.unalign}}.
}
\examples{

tau <- 10

nz <- 4
nx <- 3

Z <- matrix(1:tau, tau, nz)
X <- matrix(1:tau, tau, nx)
lags <- c(-1)
Z.na <- matrix(FALSE, tau, nz)
tt.rng <- 3:tau

xap <- als.prepare(X, Z, lags, Z.na, tt.rng, rgr = NULL, assignToEnvironment = NULL)
xap

Z[xap$y.sub.ndx, ]


############ another example

tau <- 10

nz <- 4
nx <- 3

Z <- matrix(1:tau, tau, nz)
X <- matrix(1:tau, tau, nx)
lags <- c(-1, 0, 1)
Z.na <- matrix(FALSE, tau, nz)
tt.rng <- 3:tau

xap <- als.prepare(X, Z, lags, Z.na, tt.rng, rgr = NULL, assignToEnvironment = NULL)
xap

Z[xap$y.sub.ndx, ]




## The function is currently defined as
function (X, Z, lags, Z.na, tt.rng, rgr = NULL, assignToEnvironment = NULL) 
{
    yqq <- (tt.rng[1] + min(lags, 0))
    y.tt.rng <- I(yqq:(yqq + length(tt.rng) - 1))
    if (is.null(X)) {
        xalgn <- ar.align(Z, lags, Z.na)
    }
    else {
        xalgn <- xy.align(X, Z, lags, Z.na)
    }
    X.sub <- cbind(xalgn$X.sub, rgr[xalgn$y.sub.ndx, , drop = FALSE])
    if (is.null(assignToEnvironment)) {
        return(list(X.sub = X.sub, Y.sub = xalgn$Y.sub, Y.na = xalgn$Y.na, 
            y.tt.rng = y.tt.rng, y.sub.ndx = xalgn$y.sub.ndx))
    }
    else {
        assign("X.sub", X.sub, pos = assignToEnvironment)
        assign("Y.sub", xalgn$Y.sub, pos = assignToEnvironment)
        assign("Y.na", xalgn$Y.na, pos = assignToEnvironment)
        assign("y.tt.rng", y.tt.rng, pos = assignToEnvironment)
        assign("y.sub.ndx", xalgn$y.sub.ndx, pos = assignToEnvironment)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
