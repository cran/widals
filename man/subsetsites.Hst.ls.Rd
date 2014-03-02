\name{subsetsites.Hst.ls}
\alias{subsetsites.Hst.ls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Site-Wise Extract Space-Time Covariates
}
\description{
Extract space-time covariates by site
}
\usage{
subsetsites.Hst.ls(Hst.ls, xmask)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hst.ls}{
Space-time covariates.  A list of length \eqn{\tau}, each element should be a \eqn{n} x \eqn{p_st} numeric matrix. 
}
  \item{xmask}{
Which sites to remove from \code{Hst.ls}.  A boolean vector of length \eqn{n}, or a vector of spacial indices.
}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Space-time covariates.  A list of length \eqn{\tau}, each element a \eqn{c} x \eqn{p_st} numeric matrix, where \eqn{c} is the number of \code{TRUE}'s in boolean \code{xmask}, or length of index \code{xmask}.
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
	
tau <- 70
n <- 28

Hst.ls <- list()
for(i in 1:tau) { Hst.ls[[i]] <- matrix(rnorm(n*4), nrow=n) }

subsetsites.Hst.ls(Hst.ls, c(1,3,10))


subsetsites.Hst.ls(Hst.ls, c(TRUE, TRUE, rep(FALSE, n-2)))


## The function is currently defined as
function (Hst.ls, xmask) 
{
    tau <- length(Hst.ls)
    Hst.ls.out <- list()
    for (i in 1:tau) {
        Hst.ls.out[[i]] <- Hst.ls[[i]][xmask, , drop = FALSE]
    }
    return(Hst.ls.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
