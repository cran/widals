\name{rm.cols.Hst.ls}
\alias{rm.cols.Hst.ls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove Space-Time Covariates from Model
}
\description{
Remove spacial covariates from space-time covariate list
}
\usage{
rm.cols.Hst.ls(Hst.ls, rm.col.ndx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hst.ls}{
Space-time covariates (of supporting sites).  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{pst} matrix. 
}
  \item{rm.col.ndx}{
Which columns of \code{Hst.ls} to remove.  A positive scalar integer.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
An unnamed list of length \eqn{\tau}, each element will be a numeric \eqn{n} x \eqn{pst - prm} matrix, where \eqn{prm} is the length of \code{rm.col.ndx}.
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
	
tau <- 21
n <- 7
	
pst <- 5
Hst.ls <- list()
for(i in 1:tau) { Hst.ls[[i]] <- matrix(1:pst, n, pst, byrow=TRUE) }

rm.cols.Hst.ls(Hst.ls, c(1,3))


## The function is currently defined as
function (Hst.ls, rm.col.ndx) 
{
    tau <- length(Hst.ls)
    for (i in 1:tau) {
        Hst.ls[[i]] <- Hst.ls[[i]][, -rm.col.ndx, drop = FALSE]
    }
    return(Hst.ls)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
