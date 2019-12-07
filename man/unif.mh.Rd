\name{unif.mh}
\alias{unif.mh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local Search Function
}
\description{
Search function
}
\usage{
unif.mh(n, center, sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Sample size.  A positive scalar integer.
}
  \item{center}{
Mean.  A numeric scalar (or vector).
}
  \item{sd}{
Standard deviation.  A numeric scalar (or vector).
}
}


%\details{
%%%%
%}

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
\code{\link{dlog.norm}}, \code{\link{MSS.snow}}.
}
\examples{

x <- unif.mh(100, 1, 1)
hist(x)

## The function is currently defined as
function (n, center, sd) 
{
    w <- sd * sqrt(3)
    a <- center - w
    b <- center + w
    x <- runif(n, a, b)
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
