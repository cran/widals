\name{xy.align}
\alias{xy.align}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align Data
}
\description{
Flatten, or reflatten data for lagged auto-regressive ALS
}
\usage{
xy.align(X, Y, lags, Y.mask = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Covariates.
}
  \item{Y}{
Responses.
}
  \item{lags}{
Temporal lags.
}
  \item{Y.mask}{
Indicator of missing values in \code{Y}.
}
}
\details{
This function is somewhat vestigial.  It has no meaningful role in WIDALS, but can be helpful aligning data for the sake of auto-regressive ALS.
}
\value{
A named list.
\item{X.sub}{Covariates aligned as prescribed by \code{lags}.}
\item{Y.sub}{Responses aligned as prescribed by \code{lags}.}
\item{Y.na}{\code{Y.mask} aligned as prescribed by \code{lags}.}
\item{y.sub.ndx}{The indices that map the rows of \code{Y} to the rows of \code{Y.sub}}
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
\code{\link{ar.align}}, \code{\link{als.prepare}}, \code{\link{y.unalign}}.
}
\examples{

tau <- 10

nz <- 4
nx <- 3

Y <- matrix(1:tau, tau, nz)
X <- matrix(1:tau, tau, nx)
lags <- c(-1)
Y.mask <- matrix(FALSE, tau, nz)
tt.rng <- 3:tau

xya <- xy.align(X, Y, lags, Y.mask)
xya

Y[xya$y.sub.ndx, ]


############ another example

tau <- 10

nz <- 4
nx <- 3

Y <- matrix(1:tau, tau, nz)
X <- matrix(1:tau, tau, nx)
lags <- c(-1, 0, 1)
Y.mask <- matrix(FALSE, tau, nz)
tt.rng <- 3:tau

xya <- xy.align(X, Y, lags, Y.mask)
xya

Y[xya$y.sub.ndx, ]



## The function is currently defined as
function (X, Y, lags, Y.mask = NULL) 
{
    tau <- nrow(Y)
    n <- ncol(X)
    low.ndx <- min(min(lags), 0)
    top.ndx <- max(0, max(lags))
    tau0 <- tau - (top.ndx - low.ndx)
    y.sub.ndx <- (1 - low.ndx):(tau0 - low.ndx)
    Y.sub <- Y[y.sub.ndx, , drop = FALSE]
    if (is.null(Y.mask)) {
        Y.na <- NULL
    }
    else {
        Y.na <- Y.mask[y.sub.ndx, , drop = FALSE]
    }
    X.sub <- matrix(0, tau0, n * length(lags))
    for (i in 1:length(lags)) {
        this.lag <- lags[i]
        x.start <- this.lag - low.ndx + 1
        col.rng <- ((i - 1) * n + 1):(i * n)
        X.sub[, col.rng] <- X[x.start:(tau0 + x.start - 1), , 
            drop = FALSE]
    }
    return(list(X.sub = X.sub, Y.sub = Y.sub, Y.na = Y.na, y.sub.ndx = y.sub.ndx))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
