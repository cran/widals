\name{distance}
\alias{distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spacial Distance
}
\description{
Calculate spacial distance between two sets of locations (in two-space)
}
\usage{
distance(locs1, locs2, geodesic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locs1}{
First set of locations.  E.g., supporting sites: An \eqn{n} x \code{2} matrix.  If \code{geodesic} is set to true, make sure to place latitude in first column.
}
  \item{locs2}{
Second set of locations.  E.g., interpolation sites: An \eqn{n}* x \code{2} matrix.  If \code{geodesic} is set to true, make sure to place latitude in first column.
}
  \item{geodesic}{
Use geodesic distance?  Boolean.
}
}
\details{
If \code{geodesic} is set to \code{FALSE}, Euclidean distance is returned; if \code{TRUE}, Earth's geodesic distance is returned in units kilometers.
}
\value{
An \eqn{n} x \eqn{n}* matrix.
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

locs1 <- cbind( c(-1, -1, 1, 1), c(-1, 1, -1, 1) )
locs2 <- cbind( c(0), c(0) )

distance(locs1, locs2)


locs1 <- cbind( c(32, 0), c(-114, -114) )
locs2 <- cbind( c(0), c(0) )

distance(locs1, locs2, TRUE)


####### separation of one deg long at 88 degs lat (near North-Pole) is (appx)
locs1 <- cbind( c(88), c(-114) )
locs2 <- cbind( c(88), c(-115) )
distance(locs1, locs2, TRUE)

####### separation of one deg long at 0 degs lat (Equator) is (appx)
locs1 <- cbind( c(0), c(-114) )
locs2 <- cbind( c(0), c(-115) )
distance(locs1, locs2, TRUE)





## The function is currently defined as
function (locs1, locs2, geodesic = FALSE) 
{
#    dyn.load("~/Files/Creations/C/distance.so")
    n1 <- nrow(locs1)
    n2 <- nrow(locs2)
    d.out <- rep(0, n1 * n2)
    if (geodesic) {
        D.Mx <- .C("distance_geodesic_AB", as.double(locs1[, 
            1] * pi/180), as.double(locs1[, 2] * pi/180), as.double(locs2[, 
            1] * pi/180), as.double(locs2[, 2] * pi/180), as.double(d.out), 
            as.integer(n1), as.integer(n2))[[5]]
    }
    else {
        D.Mx <- .C("distance_AB", as.double(locs1[, 1]), as.double(locs1[, 
            2]), as.double(locs2[, 1]), as.double(locs2[, 2]), 
            as.double(d.out), as.integer(n1), as.integer(n2))[[5]]
    }
    D.out <- matrix(D.Mx, n1, n2)
    return(D.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
