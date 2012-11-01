\name{load.Hst.ls.Z}
\alias{load.Hst.ls.Z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load Observations into Space-Time Covariates
}
\description{
Insert an observation matrix into space-time covariates
}
\usage{
load.Hst.ls.Z(Z, Hst.ls.Z, xwhich, rgr.lags = c(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Observation data.  A \eqn{\tau} x \eqn{n} numeric matrix.
}
  \item{Hst.ls.Z}{
Space-time covariates.  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{pst} matrix. 
}
  \item{xwhich}{
Which column of \code{Hst.ls.Z[[i]]} to insert into the \eqn{i}th row of \code{Z}.  A scalar positive integer.
}
  \item{rgr.lags}{
Temporal lagging of \code{Z}.  A scalar integer.
}
}

\details{
This function, along with \code{\link{load.Hst.ls.2Zs}}, allows the user to convert a set of observations into covariates for another set of observations.  
}

\value{
An unnamed list of length \eqn{\tau}, each element will be a numeric \eqn{n} x \eqn{pst} matrix.
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
\code{\link{load.Hst.ls.2Zs}}.
}
\examples{

###### here's an itty-bitty example

tau <- 7
n <- 5

Z <- matrix(1, tau, n)

Hst.ls <- list()
for(i in 1:tau) { Hst.ls[[i]] <- matrix(rnorm(n*4), nrow=n) }

load.Hst.ls.Z(Z, Hst.ls.Z=Hst.ls, 1, 0)


########## insert into col 3

load.Hst.ls.Z(Z, Hst.ls.Z=Hst.ls, 3, 0)




############ lag Z examples

Z <- matrix(1:tau, tau, n)

######### lag -1 Z

load.Hst.ls.Z(Z, Hst.ls.Z=Hst.ls, 1, -1)

######### lag 0 Z -- default

load.Hst.ls.Z(Z, Hst.ls.Z=Hst.ls, 1, 0)

######### lag +1 Z

load.Hst.ls.Z(Z, Hst.ls.Z=Hst.ls, 1, +1)





## The function is currently defined as
function (Z, Hst.ls.Z, xwhich, rgr.lags = c(0)) 
{
    tau <- nrow(Z)
    min.ndx <- max(1, -min(rgr.lags) + 1)
    max.ndx <- min(tau, tau - max(rgr.lags))
    for (i in min.ndx:max.ndx) {
        Hst.ls.Z[[i]][, xwhich] <- t(Z[i + rgr.lags, ])
    }
    return(Hst.ls.Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
