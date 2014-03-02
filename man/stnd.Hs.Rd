\name{stnd.Hs}
\alias{stnd.Hs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize Spacial Covariates
}
\description{
Standardize spacial covariates with respect to both the space and time dimensions
}
\usage{
stnd.Hs(Hs, Hs0 = NULL, intercept = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hs}{
Spacial covariates (of supporting sites).  An \eqn{n} x \eqn{p_s} numeric matrix.
}
  \item{Hs0}{
Spacial covariates (of interpolation sites).  An \eqn{n}* x \eqn{p_s} numeric matrix.
}
  \item{intercept}{
Include intercept term?  Boolean.
}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A named list.
\item{sHs}{An \eqn{n} x \eqn{p_s} numeric matrix.}
\item{sHs0}{An \eqn{n}* x \eqn{p_s} numeric matrix.}
\item{h.mean}{The covariates' mean over space.}
\item{h.sd}{The covariates' standard deviation over space.}
\item{n}{Number of support sites.}
\item{intercept}{The supplied intercept argument.}
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
\code{\link{stnd.Ht}}, \code{\link{stnd.Hst.ls}}, \code{\link{applystnd.Hs}}.
}
\examples{

##### Please see the examples in Hst.sumup



## The function is currently defined as
function (Hs, Hs0 = NULL, intercept = TRUE) 
{
    n <- nrow(Hs)
    h.mean <- apply(Hs, 2, mean)
    h.sd <- apply(t(t(Hs) - h.mean), 2, function(x) {
        sqrt(sum(x^2))
    })
    h.sd[h.sd == 0] <- 1
    sHs <- t((t(Hs) - h.mean)/h.sd)
    if (intercept) {
        sHs[, 1] <- 1/sqrt(n)
    }
    sHs0 <- NULL
    if (!is.null(Hs0)) {
        sHs0 <- t((t(Hs0) - h.mean)/h.sd)
        if (intercept) {
            sHs0[, 1] <- 1/sqrt(n)
        }
    }
    ls.out <- list(sHs = sHs, sHs0 = sHs0, h.mean = h.mean, h.sd = h.sd, 
        n = n, intercept = intercept)
    return(ls.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
