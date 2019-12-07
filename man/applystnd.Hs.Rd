\name{applystnd.Hs}
\alias{applystnd.Hs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize Spacial Covariates with Existing Object 
}
\description{
Standardize spacial covariates with respect to both the space and time dimensions
}
\usage{
applystnd.Hs(Hs0, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hs0}{
Spacial covariates (of interpolation sites).  An \eqn{n}* x \eqn{p_s} numeric matrix.
}
  \item{x}{
Spacial standardization object, as created by \code{\link{stnd.Hs}}.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
An \eqn{n}* x \eqn{p_s} matrix.
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
\code{\link{stnd.Hst.ls}}, \code{\link{applystnd.Hst.ls}}.
}


\examples{

n.all <- 21
Hs.all <- cbind(1, rnorm(n.all, 1, 0.1), rnorm(n.all, -200, 21))

ndx.interp <- c(1,3,5)
ndx.support <- I(1:n.all)[ -ndx.interp ]


Hs <- Hs.all[ndx.support, , drop=FALSE]

xsns.obj <- stnd.Hs(Hs)

Hs0 <- Hs.all[ndx.interp, , drop=FALSE]

sHs0 <- applystnd.Hs(Hs0, xsns.obj)
sHs0

xsns.obj$sHs

crossprod(xsns.obj$sHs) / nrow(Hs)

crossprod(sHs0) / nrow(sHs0)


## The function is currently defined as
function (Hs0, x) 
{
    sHs0 <- t((t(Hs0) - x$h.mean)/x$h.sd)
    if (x$intercept) {
        sHs0[, 1] <- 1/sqrt(x$n)
    }
    return(sHs0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
