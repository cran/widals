\name{dlog.norm}
\alias{dlog.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local Search Function
}
\description{
Local hyperparameter exponentiated-normal search function
}
\usage{
dlog.norm(n, center, sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size.  A positive scalar integer.
}
  \item{center}{
Exponential of the mean.  A numeric scalar (or vector).
}
  \item{sd}{
Standard deviation.  A numeric scalar (or vector).
}
}
\details{
This function can be used by \code{\link{MSS.snow}}.
}
\value{
A numeric vector of length \eqn{n}.
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
\code{\link{unif.mh}}, \code{\link{MSS.snow}}, \code{\link{fun.load}}.
}

\examples{

x <- dlog.norm(100, 1, 1)
hist(x)

## The function is currently defined as
function (n, center, sd) 
{
    return(exp(rnorm(n, log(center), sd)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
