\name{stnd.Hst.ls}
\alias{stnd.Hst.ls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize Space-Time Covariates
}
\description{
Standardize spacio-temporal covariates with respect to both the spacial and time dimensions
}
\usage{
stnd.Hst.ls(Hst.ls, Hst0.ls = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hst.ls}{
Space-time covariates (of supporting sites).  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{p_st} matrix. 
}
  \item{Hst0.ls}{
Space-time covariates (of interpolation sites).  A list of length \eqn{\tau}, each element should be a numeric \eqn{n}* x \eqn{p_st} matrix. 
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A named list.
\item{sHst.ls}{A list of length \eqn{\tau}, each element a numeric \eqn{n} x \eqn{p_st} matrix.}
\item{sHst0.ls}{A list of length \eqn{\tau}, each element a \eqn{n}* x \eqn{p_st} matrix}
\item{h.mean}{The covariates' mean over space-time.}
\item{h.sd}{The covariates' standard deviation over space-time.}
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
\code{\link{stnd.Ht}}, \code{\link{stnd.Hs}}, \code{\link{applystnd.Hst.ls}}.
}
\examples{
	
##### Please see the examples in Hst.sumup

## The function is currently defined as
function (Hst.ls, Hst0.ls = NULL) 
{
    tau <- length(Hst.ls)
    big.sum <- 0
    for (i in 1:tau) {
        big.sum <- big.sum + apply(Hst.ls[[i]], 2, mean)
    }
    h.mean <- big.sum/tau
    sHst.ls <- list()
    big.sum.mx <- 0
    for (i in 1:tau) {
        sHst.ls[[i]] <- t(t(Hst.ls[[i]]) - h.mean)
        big.sum.mx <- big.sum.mx + crossprod(sHst.ls[[i]])
    }
    cov.mx <- big.sum.mx/tau
    sqrtXX <- 1/sqrt(diag(cov.mx))
    for (i in 1:tau) {
        sHst.ls[[i]] <- t(t(sHst.ls[[i]]) * sqrtXX)
    }
    sHst0.ls <- NULL
    if (!is.null(Hst0.ls)) {
        sHst0.ls <- list()
        for (i in 1:tau) {
            sHst0.ls[[i]] <- t((t(Hst0.ls[[i]]) - h.mean) * sqrtXX)
        }
    }
    ls.out <- list(sHst.ls = sHst.ls, sHst0.ls = sHst0.ls, h.mean = h.mean, 
        h.sd = 1/sqrtXX)
    return(ls.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
