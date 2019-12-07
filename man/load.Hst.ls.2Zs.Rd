\name{load.Hst.ls.2Zs}
\alias{load.Hst.ls.2Zs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load Observations into Space-Time Covariates
}
\description{
Insert an observation matrix into space-time covariates, but segregate based on missing values
}
\usage{
load.Hst.ls.2Zs(Z, Z.na, Hst.ls.Z, xwhich, rgr.lags = c(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Observation data.  A \eqn{\tau} x \eqn{n} numeric matrix.
}
  \item{Z.na}{
Missing data indicator.  A \eqn{\tau} x \eqn{n} boolean matrix.
}
  \item{Hst.ls.Z}{
Space-time covariates.  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{p_st} matrix.
}
  \item{xwhich}{
Which column-pair of \code{Hst.ls.Z[[i]]} to insert into the \eqn{i}th row of \code{Z}.  A scalar positive integer.  By 'column-pair', we mean, e.g., a value of 1 will fill columns 1 and 2, a value of 2 will fill columns 3 and 4, a value of 3 will fill columns 5 and 6, etc.
}
  \item{rgr.lags}{
Temporal lagging of \code{Z}.  A scalar integer.
}

}
\details{
This function, along with \code{\link{load.Hst.ls.Z}}, allows the user to convert a set of observations into covariates for another set of observations.  Unlike \code{\link{load.Hst.ls.Z}}, this function \emph{splits} \code{Z} based on the argument \code{Z.na}.  Values associated with \code{FALSE} elements of \code{Z.na} are placed into the first column of the specified column-pair of \code{Hst.ls.Z}, Values associated with \code{TRUE} elements of \code{Z.na} are placed into the second column of the specified column-pair of \code{Hst.ls.Z} (all other values in in the specified column-pair of \code{Hst.ls.Z} are zeroed).
}
\value{
An unnamed list of length \eqn{\tau}, each element will be a numeric \eqn{n} x \eqn{p_st} matrix.
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
\code{\link{load.Hst.ls.Z}}.
}


\examples{


###### here's an itty-bitty example

tau <- 7
n <- 5

Z <- matrix(1, tau, n)

Z.na <- matrix(FALSE, tau, n)
Z.na[2:3, 4] <- TRUE

Z[Z.na] <- 2

Hst.ls <- list()
for(i in 1:tau) { Hst.ls[[i]] <- matrix(rnorm(n*4), nrow=n) }


load.Hst.ls.2Zs(Z, Z.na, Hst.ls.Z=Hst.ls, 1, 0)


########## insert into cols 3 and 4

load.Hst.ls.2Zs(Z, Z.na, Hst.ls.Z=Hst.ls, 2, 0)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
