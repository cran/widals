\name{Hst.sumup}
\alias{Hst.sumup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Covariance Matrix
}
\description{
Calculate the covariance matrix of all model covariates
}
\usage{
Hst.sumup(Hst.ls, Hs = NULL, Ht = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hst.ls}{
Space-time covariates.  A list of length \eqn{\tau}, each element containing a \eqn{n} x \eqn{pst} numeric matrix.
}
  \item{Hs}{
Spacial covariates.  An \eqn{n} x \eqn{ps} numeric matrix.
}
  \item{Ht}{
Temporal covariates.  An \eqn{\tau} x \eqn{pt} numeric matrix.
}
}
\details{
Important: The order of the arguments in this function is NOT the same as in the returned covariance matrix.  The order in the covariance matrix is the same as in other functions in this package:  \code{Hs}, \code{Ht}, \code{Hst.ls}.
}
\value{
A \eqn{(ps+pt+pst)} x \eqn{(ps+pt+pst)} numeric, symmetrix, non-negative definite matrix.
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
	
tau <- 20
n <- 10
Ht <- cbind(sin(1:tau), cos(1:tau))

Hs <- cbind(rnorm(10), rnorm(n, 5, 49))

Hst.ls <- list()
for(tt in 1:tau) {
Hst.ls[[tt]] <- cbind(rnorm(n, 1, 0.1), rnorm(n, -200, 21))
}


Hst.sumup(Hst.ls, Hs, Ht)



########### standardize all covariates

x1 <- stnd.Hst.ls(Hst.ls, NULL)$sHst.ls
x2 <- stnd.Hs(Hs, NULL, FALSE)$sHs
x3 <- stnd.Ht(Ht, n)


Hst.sumup(x1, x2, x3)



## The function is currently defined as
function (Hst.ls, Hs = NULL, Ht = NULL) 
{
    tau <- length(Hst.ls)
    if(tau < 1) { tau <- nrow(Ht) }
    if(is.null(tau)) { tau <- 10 ; cat("tau assumed to be 10.", "\n") }
    n <- nrow(Hst.ls[[1]])
    if(is.null(n)) { n <- nrow(Hs) }
    big.sum <- 0
    for (i in 1:tau) {
        if (!is.null(Ht)) {
            Ht.mx <- matrix(Ht[i, ], n, ncol(Ht), byrow = TRUE)
        }
        else {
            Ht.mx <- NULL
        }
        big.sum <- big.sum + crossprod(cbind(Hs, Ht.mx, Hst.ls[[i]]))
    }
    return(big.sum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
