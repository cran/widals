\name{H.Earth.solar}
\alias{H.Earth.solar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solar Radiation
}
\description{
Calculate Incident Solar Area (ISA)
}
\usage{
H.Earth.solar(x, y, dateDate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Longitude.  Numeric vector of length \eqn{n}.
}
  \item{y}{
Latitude.  Numeric vector of length \eqn{n}.
}
  \item{dateDate}{
Posix date.  Numeric vector of length \eqn{\tau}.
}
}
\details{
This function returns a spacio-temporal covariate list (Earth's ISA is space-time \emph{non-seperable}).  A negative value indicates that at that time (list index), and at that location (matrix row), the sun is below the horizon all day.
}
\value{
An unnamed list of length \eqn{\tau}, each element of which is an \eqn{n} x 1 matrix.  
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

lat <- c(0, -88)
lon <- c(0, 0)
dateDate <- strptime( c('20120621', '20120320'), '\%Y\%m\%d')

H.Earth.solar(lon, lat, dateDate)


## The function is currently defined as
function (x, y, dateDate) 
{
    Hst.ls <- list()
    n <- length(y)
    tau <- length(dateDate)
    equinox <- strptime("20110320", "\%Y\%m\%d")
    for (i in 1:tau) {
        this.date <- dateDate[i]
        dfe <- as.integer(difftime(this.date, equinox, units = "day"))
        dfe
        psi <- 23.5 * sin(2 * pi * dfe/365.25)
        psi
        eta <- 90 - (360/(2 * pi)) * acos(cos(2 * pi * y/360) * 
            cos(2 * pi * psi/360) + sin(2 * pi * y/360) * sin(2 * 
            pi * psi/360))
        surface.area <- sin(2 * pi * eta/360)
        surface.area
        Hst.ls[[i]] <- cbind(surface.area)
    }
    return(Hst.ls)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
