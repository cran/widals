\name{create.rm.ndx.ls}
\alias{create.rm.ndx.ls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-Validation Indices
}
\description{
Create a list of vectors of indices to remove for \emph{k}-fold cross-validation
}
\usage{
create.rm.ndx.ls(n, xincmnt = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of sites.  A scalar integer.
}
  \item{xincmnt}{
How many cv folds, i.e., \emph{k}.
}
}
\details{
The name of the object produced by this function is commonly \code{rm.ndx} in this documentation.  See \code{\link{MSS.snow}} for a reminder that this object is passed out-of-scope when using \code{MSS.snow}.

In this package \code{rm.ndx} is used by \code{\link{Hals.fastcv.snow}} and \code{\link{widals.snow}}; however, creating this object as a list using this function is only necessary when using \code{\link{widals.snow}} with \code{cv=2} (i.e., 'true' cross-validation).
}
\value{
An unnamed list of integer (>0) vectors.
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
n <- 100
xincmnt <- 7
rm.ndx <- create.rm.ndx.ls(n=n, xincmnt=xincmnt)
rm.ndx

######## if we want randomization of indices:
n <- 100
xincmnt <- 7
rm.ndx <- create.rm.ndx.ls(n=n, xincmnt=xincmnt)

rnd.ndx <- sample(I(1:n))
for(i in 1:length(rm.ndx)) { rm.ndx[[i]] <- rnd.ndx[rm.ndx[[i]]] }
rm.ndx

## The function is currently defined as
function (n, xincmnt = 10) 
{
    rm.ndx.ls <- list()
    for (i in 1:xincmnt) {
        xrm.ndxs <- seq(i, n + xincmnt, by = xincmnt)
        xrm.ndxs <- xrm.ndxs[xrm.ndxs <= n]
        rm.ndx.ls[[i]] <- xrm.ndxs
    }
    return(rm.ndx.ls)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
