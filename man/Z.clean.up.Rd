\name{Z.clean.up}
\alias{Z.clean.up}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clean Data
}
\description{
A crude, brute-force way to destroy bad values in data.
}
\usage{
Z.clean.up(Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Data.  A \eqn{\tau} x \eqn{n} matrix.
}
}
\details{
This function replaces intractable values, e.g., \code{NA}, or \code{-Inf}, in data, with the global mean. 
}
\value{
A \eqn{\tau} x \eqn{n} numeric matrix.
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

tau <- 10
n <- 7

Z <- matrix(1, tau, n)
Z[2,4] <- -Inf
Z[3,4] <- Inf
Z[4,4] <- NA
Z[5,4] <- log(-1)
Z

Z.clean.up(Z)





## The function is currently defined as
function (Z) 
{
    Z[Z == Inf | Z == -Inf] <- NA
    Z[is.na(Z) | is.nan(Z)] <- mean(Z, na.rm = TRUE)
    return(Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
