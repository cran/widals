\name{O3}
\alias{O3}
\docType{data}
\title{
California Ozone
}
\description{
Daily airborne Ozone concentrations (ppb) over California, 68 fixed sensors, 2005-2006
}
\usage{data(O3)}
\format{
  The format is:
  
List of 5

\code{$Z}      :'data.frame': 730 obs. of  68 variables: Ozone ppb for 68 sensors.
 
\code{$locs}   :'data.frame': 68 obs. of  2 variables: longitude, latitude for 68 sites.
 
\code{$helevs} : elevations (meters) for 68 sites.
 
\code{$locs0}  : 2358 obs. of  2 variables: longitude, latitude for interpolation grid.
 
\code{$helevs0}: elevations (meters) for 2358 interpolation grid sites.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
The Ozone data originates from the California Air Resources Board (CARB).  The interpolation grid elevations originate from the Google Elevation API.
}
\references{
The O3 data comes from the California Air Resources Board (CARB).
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(O3)
}
\keyword{datasets}
