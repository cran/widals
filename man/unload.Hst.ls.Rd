\name{unload.Hst.ls}
\alias{unload.Hst.ls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a Space-Time Covariate into Data
}
\description{
Convert a spacio-temporal covariate into contemporaneous data
}
\usage{
unload.Hst.ls(Hst.ls, which.col, rgr.lags)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hst.ls}{
Space-time covariates.  A list of length \eqn{\tau}, each element should be a \eqn{n} x \eqn{p_st} numeric matrix.
}
  \item{which.col}{
Which column of \code{Hst.ls[[i]]} to insert into the \eqn{i}th row of \code{Z}.  A scalar positive integer. 
}
  \item{rgr.lags}{
Temporal lagging of \code{Z}.  A scalar integer.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A numeric \eqn{\tau} x \eqn{n} matrix.
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
\code{\link{load.Hst.ls.Z}}, \code{\link{load.Hst.ls.2Zs}}.
}
\examples{
	
###### here's an itty-bitty example

tau <- 7
n <- 5

Hst.ls <- list()
for(i in 1:tau) { Hst.ls[[i]] <- matrix(rnorm(n*4), nrow=n) }

Zh <- unload.Hst.ls(Hst.ls, 1, 0)


## The function is currently defined as
function (Hst.ls, which.col, rgr.lags) 
{
    n <- nrow(Hst.ls[[1]])
    tau <- length(Hst.ls)
    Z.out <- matrix(NA, tau, n)
    min.ndx <- max(1, -min(rgr.lags) + 1)
    max.ndx <- min(tau, tau - max(rgr.lags))
    for (i in min.ndx:max.ndx) {
        Z.out[i - rgr.lags, ] <- Hst.ls[[i]][, which.col]
    }
    return(Z.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
