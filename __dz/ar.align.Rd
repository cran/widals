\name{ar.align}
\alias{ar.align}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align Data
}
\description{
Flatten, or reflatten data for lagged auto-regressive ALS
}
\usage{
ar.align(Z, lags, Z.mask = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Responses.
}
  \item{lags}{
Temporal lags.
}
  \item{Z.mask}{
Indicator of missing values in \code{Z}.
}
}
\details{
This function is somewhat vestigial.  It has no meaningful role in WIDALS, but can be helpful aligning data for the sake of auto-regressive ALS.
}
\value{
A named list.
\item{X.sub}{Responses re-aligned as predictors as prescribed by \code{lags}.}
\item{Y.sub}{Responses aligned as prescribed by \code{lags}.}
\item{Y.na}{\code{Z.na} aligned as prescribed by \code{lags}.}
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
\code{\link{als.prepare}}, \code{\link{xy.align}}, \code{\link{y.unalign}}.
}
\examples{

tau <- 10

nz <- 4
nx <- 3

Z <- matrix(1:tau, tau, nz)
lags <- c(-1)
Z.na <- matrix(FALSE, tau, nz)

xara <- ar.align(Z, lags, Z.na)
xara

Z[xara$y.sub.ndx, ]


############ another example

tau <- 10

nz <- 4
nx <- 3

Z <- matrix(1:tau, tau, nz)
lags <- c(-1, 0, 1)
Z.na <- matrix(FALSE, tau, nz)


xara <- ar.align(Z, lags, Z.na)
xara

Z[xara$y.sub.ndx, ]



## The function is currently defined as
function (Z, lags, Z.mask = NULL) 
{
    tau <- nrow(Z)
    n <- ncol(Z)
    low.ndx <- min(min(lags), 0)
    top.ndx <- max(0, max(lags))
    tau0 <- tau - (top.ndx - low.ndx)
    y.sub.ndx <- (1 - low.ndx):(tau0 - low.ndx)
    Y <- Z[y.sub.ndx, , drop = FALSE]
    if (is.null(Z.mask)) {
        Y.na <- NULL
    }
    else {
        Y.na <- Z.mask[y.sub.ndx, , drop = FALSE]
    }
    X <- matrix(0, tau0, n * length(lags))
    for (i in 1:length(lags)) {
        this.lag <- lags[i]
        x.start <- this.lag - low.ndx + 1
        col.rng <- ((i - 1) * n + 1):(i * n)
        X[, col.rng] <- Z[x.start:(tau0 + x.start - 1), , drop = FALSE]
    }
    return(list(X.sub = X, Y.sub = Y, Y.na = Y.na, y.sub.ndx = y.sub.ndx))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
