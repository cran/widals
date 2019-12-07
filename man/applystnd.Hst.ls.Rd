\name{applystnd.Hst.ls}
\alias{applystnd.Hst.ls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize Space-Time Covariates with Existing Object 
}
\description{
Standardize spacio-temporal covariates with respect to both the space and time dimensions
}
\usage{
applystnd.Hst.ls(Hst0.ls, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hst0.ls}{
Space-time covariates (of interpolation sites).  A list of length \eqn{\tau}, each element should be a \eqn{n}* x \eqn{p_st} numeric matrix.
}
  \item{x}{
Space-time standardization object, as created by \code{\link{stnd.Hst.ls}}.
}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
	An unnamed list of length \eqn{\tau}, each element a \eqn{n}* x \eqn{p_st} numeric matrix. 
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
n.all <- 10
	
Hst.ls.all <- list()
for(tt in 1:tau) {
	Hst.ls.all[[tt]] <- cbind(rnorm(n.all, 1, 0.1), rnorm(n.all, -200, 21))
}

ndx.interp <- c(1,3,5)
ndx.support <- I(1:n.all)[ -ndx.interp ]

Hst.ls <- subsetsites.Hst.ls(Hst.ls.all, ndx.support)

xsnst.obj <- stnd.Hst.ls(Hst.ls)

Hst0.ls <- subsetsites.Hst.ls(Hst.ls.all, ndx.interp)

sHst0.ls <- applystnd.Hst.ls(Hst0.ls, xsnst.obj)



Hst.sumup(xsnst.obj$sHst.ls)

Hst.sumup(sHst0.ls)



## The function is currently defined as
function (Hst0.ls, x) 
{
    tau <- length(Hst0.ls)
    sHst0.ls <- list()
    for (i in 1:tau) {
        sHst0.ls[[i]] <- t((t(Hst0.ls[[i]]) - x$h.mean)/x$h.sd)
    }
    return(sHst0.ls)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
