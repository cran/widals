\name{stnd.Ht}
\alias{stnd.Ht}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize Temporal Covariates
}
\description{
Standardize temporal covariates with respect to both the spacial and time dimensions
}
\usage{
stnd.Ht(Ht, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ht}{
Temporal covariates (of supporting sites).  A \eqn{\tau} x \eqn{pt} numeric matrix.
}
  \item{n}{
Number of sites.  A positive scalar integer.
}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \eqn{\tau} x \eqn{pt} numeric matrix.
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
\code{\link{stnd.Hs}}, \code{\link{stnd.Hst.ls}}.
}
\examples{

##### Please see the examples in Hst.sumup


## The function is currently defined as
function (Ht, n) 
{
    h.mean <- apply(Ht, 2, mean)
    sHt <- t(t(Ht) - h.mean)
    sHt <- t(t(sHt)/apply(sHt, 2, function(x) {
        sqrt(sum(x^2))
    }))
    sHt <- sHt * sqrt(nrow(Ht)/n)
    return(sHt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
