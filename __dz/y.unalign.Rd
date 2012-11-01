\name{y.unalign}
\alias{y.unalign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align Data
}
\description{
Flatten, or reflatten data for lagged auto-regressive ALS
}
\usage{
y.unalign(Y.sub, lags)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.sub}{
Aligned response, e.g., as returned by \code{\link{als.prepare}}.
}
  \item{lags}{
Temporal lags.
}
}
\details{
This function is somewhat vestigial.  It has no meaningful role in WIDALS, but can be helpful aligning data for the sake of auto-regressive ALS.
}
\value{
Responses stretched row-wise.  The inverse effect of the conversion of \code{Y} to \code{Y.sub} by \code{\link{xy.align}}.
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
\code{\link{ar.align}}, \code{\link{xy.align}}, \code{\link{als.prepare}}.
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

xya <- xy.align(X, Y, lags, Y.mask) #### align
xya

Y[xya$y.sub.ndx, ]

##### now, restore Y
y.unalign(xya$Y.sub, lags)




## The function is currently defined as
function (Y.sub, lags) 
{
    tau0 <- nrow(Y.sub)
    low.ndx <- min(min(lags), 0)
    top.ndx <- max(0, max(lags))
    tau <- tau0 + (top.ndx - low.ndx)
    Y <- matrix(NA, tau, ncol(Y.sub))
    Y[(1 - low.ndx):(tau0 - low.ndx), ] <- Y.sub
    return(Y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
