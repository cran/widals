\name{fuse.Hst.ls}
\alias{fuse.Hst.ls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge Contemporaneous Space-Time Covariates
}
\description{
Fuse together two lists of spacio-temporal covariates
}
\usage{
fuse.Hst.ls(Hst.ls1, Hst.ls2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hst.ls1}{
Space-time covariates.  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{p_1} matrix. 
}
  \item{Hst.ls2}{
Space-time covariates.  A list of length \eqn{\tau}, each element should be a numeric \eqn{n} x \eqn{p_2} matrix. 
}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
An unnamed list of length \eqn{\tau}, each element will be a numeric \eqn{n} x \eqn{(p_1+p_2)} matrix.
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

set.seed(9999)

tau <- 5
n <- 7

p1 <- 2
Hst.ls1 <- list()
for(i in 1:tau) { Hst.ls1[[i]] <- matrix(rnorm(n*p1), nrow=n) }

p2 <- 3
Hst.ls2 <- list()
for(i in 1:tau) { Hst.ls2[[i]] <- matrix(rnorm(n*p2), nrow=n) }

fuse.Hst.ls(Hst.ls1, Hst.ls2)


## The function is currently defined as
function (Hst.ls1, Hst.ls2) 
{
    tau <- length(Hst.ls1)
    for (i in 1:tau) {
        Hst.ls1[[i]] <- cbind(Hst.ls1[[i]], Hst.ls2[[i]])
    }
    return(Hst.ls1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
