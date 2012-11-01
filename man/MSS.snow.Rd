\name{MSS.snow}
\alias{MSS.snow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Metaheuristic Stochastic Search
}
\description{
Locate WIDALS hyperparameters
}
\usage{
MSS.snow(FUN.source, current.best, p.ndx.ls, f.d, sds.mx, k.glob, k.loc.coef, X = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FUN.source}{
Search function definitions (see Details).  A path to source code, or function, e.g., \code{\link{fun.load.widals.a}}.
}
  \item{current.best}{
An initial cost.  A scalar.  Setting to NA will cause \code{MSS.snow} to make an initial pass over the data to create an initial cost to beat.
}
  \item{p.ndx.ls}{
Hyperparameter indices (of \code{GP}) to search.  A list of vectors.  For example, \code{list( c(1,2), c(3,4,5) )} will instruct \code{MSS.snow}, for each local search, to search over the first two hyperparameters as a pair, then to search the last three as a group.
}
  \item{f.d}{
Local search functions.  A list of functions (one for each element of GP).  Typically, for WIDALS, all five will be \code{\link{dlog.norm}}.
}
  \item{sds.mx}{
The standard deviations for \code{f.d}.  An \emph{k.glob} x \emph{q} matrix, where \emph{q} is the number of hyperparameters, i.e., the length of GP.
}
  \item{k.glob}{
The number of global searches.  A scalar integer.
}
  \item{k.loc.coef}{
The coeficient for the number of local searches to make.  A scalar integer.
}
  \item{X}{
A placeholder for values to be passed between functions inside \code{MSS.snow} (see Details).
}
}
\details{
This function requires the presence of a number of values and functions out-of-scope.  It is assumed that these are available in the Global Environment.  They are: \code{run.parallel} (boolean), \code{FUN.MH} (a function that creates, for a given \code{GP}, a cost), \code{FUN.GP} (a function that applies constraints to \code{GP}), \code{FUN.I} (a function that does something when local searches have reduced the cost), \code{FUN.EXIT} (a function that does something when \code{MSS.snow} is done). 

Examine the code for \code{\link{fun.load.widals.a}} for an example of the four functions described above.  Note that these four functions may themselves require objects out-of-scope.

In general, for a given \code{R} session, special care should be taken concerning the naming and assigning of the following objects: \code{Z} (the space-time data), \code{Z.na} (a boolean matrix indicating missing values in \code{Z}), \code{locs} (site locations), \code{Hs} (spacial covariates), \code{Ht} (temporal covariates), \code{Hst.ls} (space-time covariates), \code{lags} (temporal lag vector), \code{b.lag} (the ALS lag), \code{cv} (cross-validation switch), \code{xgeodesic} (boolean), \code{ltco} (weight cut-off), \code{GP} (hyperparameter vector), \code{run.parralel} (boolean), \code{stnd.d} (boolean), \code{train.rng} (time index vector), \code{test.rng} (time index vector).
}
\value{
Nothing.  After completion, the best hyperparameters, \code{GP}, are assigned to the Global Environment.
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
\code{\link{Hals.fastcv.snow}}, \code{\link{Hals.snow}}, \code{\link{widals.snow}}.
}


\examples{

##### simulate a state-space system (using pkg SSsimple)

\dontrun{
	
tau <- 77 #### number of time points

d.alpha <- 2
R.scale <- 1
sigma2 <- 0.01
F <- 0.999
Q <- 0.1

udom <- (0:300)/100
plot( udom,    R.scale * exp(-d.alpha*udom) ,  type="l", col="red" ) #### see the covariogram

n.all <- 70 ##### number of spacial locations

set.seed(9999)
locs.all <- cbind(runif(n.all, -1, 1), runif(n.all, -1, 1)) #### random location of sensors

D.mx <- distance(locs.all, locs.all, FALSE) #### distance matrix

R.all <- exp(-d.alpha*D.mx) + diag(sigma2, n.all) #### create measurement variance using distance and covariogram

Hs.all <- matrix(1, n.all, 1) #### constant mean function

xsssim <- SS.sim(F=F, H=Hs.all, Q=Q, R=R.all, length.out=tau, beta0=0) ##### use SSsimple to simulate system

Z.all <- xsssim$Z ###### system observation matrix

	
	
########### now make assignments required by MSS.snow
	
set.seed(9999)
ndx.interp <- sample(1:n.all, size=5) ##### randomly remove five sites to serve as interpolation points
ndx.support <- I(1:n.all)[ -ndx.interp ] ##### support sites



####################################### what follows are important assignments, since MSS.snow and the four helper functions
####################################### will look for these in the Global Environment to commence fitting the model (as noted in Details above)
train.rng <- 30:(tau) ; test.rng <- train.rng

Z <- Z.all[ , ndx.support ] 
Hs <- Hs.all[ ndx.support, , drop=FALSE] 
locs <- locs.all[ndx.support, , drop=FALSE] 

Ht <- NULL
Hst.ls <- NULL

lags <- c(0) 
b.lag <- c(-1) 
cv <- -2
xgeodesic <- FALSE
stnd.d <- FALSE
ltco <- -10
GP <- c(1/10, 1, 20, 20, 1) ### -- initial hyperparameter values
run.parallel <- TRUE 

if( cv==2 ) { rm.ndx <- create.rm.ndx.ls( nrow(Hs), 14 ) } else { rm.ndx <- 1:nrow(Hs) }
rgr.lower.limit <- 10^(-7) ; d.alpha.lower.limit <- 10^(-3) ; rho.upper.limit <- 10^(4)


############## tell snowfall to use two threads for local searches
sfInit(TRUE, cpus=2)


######## now, finally, search for best fit over support -- note that p.ndx.ls and f.d are produced inside fun.load.widals.a()
MSS.snow(fun.load.widals.a, NA, p.ndx.ls, f.d, matrix(1/10, 10, length(GP)), 10, 7)
sfStop()

######## we can use these hyperparameters to interpolate to the deliberately removed sites, and measure MSE, RMSE
Z0.hat <- widals.predict(Z, Hs, Ht, Hst.ls, locs, lags, b.lag, Hs0=Hs.all[ ndx.interp, , drop=FALSE ], Hst0.ls=NULL, locs0=locs.all[ ndx.interp, , drop=FALSE],
geodesic = xgeodesic, wrap.around = NULL, GP, stnd.d = stnd.d, ltco = ltco)

resids.wid <- ( Z.all[ , ndx.interp ] - Z0.hat )
mse.wid <- mean( resids.wid[ test.rng, ]^2 )
mse.wid
sqrt(mse.wid)






########################################### Simulated Imputation with WIDALS
Z.all <- xsssim$Z
Z.missing <- Z.all

Z.na.all <- matrix( sample(c(TRUE, FALSE), size=n.all*tau, prob=c(0.01, 0.99), replace=TRUE), tau, n.all)
Z.missing[ Z.na.all ] <- NA


Z <- Z.missing
Z[ is.na(Z) ] <- mean(Z, na.rm=TRUE)
X <- list("Z.fill"=Z)

Z.na <- Z.na.all
Hs <- Hs.all
locs <- locs.all
Ht <- NULL
Hst.ls <- NULL
lags <- c(0)
b.lag <- c(-1)
cv <- -2
xgeodesic <- FALSE
ltco <- -10
if( cv==2 ) { rm.ndx <- create.rm.ndx.ls( nrow(Hs), 14 ) } else { rm.ndx <- 1:nrow(Hs) }

GP <- c(1/10, 1, 20, 20, 1)

rgr.lower.limit <- 10^(-7) ; d.alpha.lower.limit <- 10^(-3) ; rho.upper.limit <- 10^(4)

run.parallel <- TRUE

sfInit(TRUE, cpus=2)

MSS.snow(fun.load.widals.fill, NA, p.ndx.ls, f.d, seq(2, 0.01, length=10)*matrix(1/10, 10, length(GP)), 10, 7, X=X)
sfStop()

sqrt(mean(( (Z.all[train.rng, ] - Z.fill[train.rng, ])^2 )[ Z.na[ train.rng, ] ]))



    
    
    
    

############################################ Now Try with ALS alone

Z.all <- xsssim$Z

GP <- c(1/10, 1) ### -- initial hyperparameter values

############## tell snowfall to use two threads for local searches
sfInit(TRUE, cpus=2)


######## now, finally, search for best fit over support -- note that p.ndx.ls and f.d are produced inside fun.load.widals.a()
MSS.snow(fun.load.hals.a, NA, p.ndx.ls, f.d, matrix(1/10, 10, length(GP)), 10, 7)
sfStop()

######## we can use these hyperparameters to interpolate to the deliberately removed sites, and measure MSE, RMSE
hals.obj <- H.als.b(Z, Hs, Ht, Hst.ls, rho=GP[1], reg=GP[2], b.lag = b.lag, Hs0 = Hs.all[ ndx.interp, , drop=FALSE ], Ht0 = NULL, Hst0.ls = NULL)
Z0.hat <- hals.obj$Z0.hat

resids.als <- ( Z.all[ , ndx.interp ] - Z0.hat )
mse.als <- mean( resids.als[ test.rng, ]^2 )
mse.als
sqrt(mse.als)



########################################### Simulated Imputation with ALS
Z.all <- xsssim$Z
Z.missing <- Z.all

set.seed(99)
Z.na.all <- matrix( sample(c(TRUE, FALSE), size=n.all*tau, prob=c(0.03, 0.97), replace=TRUE), tau, n.all)
Z.missing[ Z.na.all ] <- NA


Z <- Z.missing
Z[ is.na(Z) ] <- 0 #mean(Z, na.rm=TRUE)
X <- list("Z.fill"=Z)
    
Z.na <- Z.na.all

Hs <- Hs.all

GP <- c(1/10, 1) ### -- initial hyperparameter values

sfInit(TRUE, cpus=2)
MSS.snow(fun.load.hals.fill, NA, p.ndx.ls, f.d, seq(3, 0.01, length=10)*matrix(1, 10, length(GP)), 10, 7, X=X)


sqrt(mean(( (Z.all[train.rng, ] - Z.fill[train.rng, ])^2 )[ Z.na[ train.rng, ] ]))

}








## The function is currently defined as
function (FUN.source, current.best, p.ndx.ls, f.d, sds.mx, k.glob, k.loc.coef, X = NULL)
{
    envmh <- environment(NULL)
    GP <- GP
    if(is.function(FUN.source)) {
        FUN.source()
    } else {
        if (!is.null(FUN.source)) {
            source(FUN.source)
        }
    }
    if (is.na(current.best)) {
        GP.mx <- matrix(GP, 1, length(GP))
        if (!is.null(FUN.GP)) {
            GP.mx <- FUN.GP(GP.mx)
        }
        current.best <- FUN.MH(1, GP.mx = GP.mx, X = X)
        cat(current.best, "\n")
    }
    if (!is.null(k.glob)) {
        for (k.times in 1:k.glob) {
            cat(k.times, "\n")
            for (p.ndx in p.ndx.ls) {
                n.mh <- as.integer(k.loc.coef * 2^length(p.ndx))
                GP.mx <- matrix(GP, n.mh, length(GP), byrow = TRUE)
                for (ip in p.ndx) {
                    GP.mx[, ip] <- f.d[[ip]](n.mh, GP[ip], sds.mx[k.times,
                    ip])
                }
                if (!is.null(FUN.GP)) {
                    GP.mx <- FUN.GP(GP.mx)
                }
                if (run.parallel) {
                    sfOut <- sfClusterApplyLB(1:n.mh, FUN.MH, GP.mx = GP.mx,
                    X = X)
                }
                else {
                    sfOut <- list()
                    for (jj in 1:n.mh) {
                        sfOut[[jj]] <- FUN.MH(jj, GP.mx = GP.mx,
                        X = X)
                    }
                }
                errs <- unlist(sfOut)
                errs[is.na(errs)] <- Inf
                errs[is.nan(errs)] <- Inf
                best.ndx <- which(errs == min(errs))[1]
                if (errs[best.ndx] < current.best) {
                    current.best <- errs[best.ndx]
                    GP <- GP.mx[best.ndx, , drop = TRUE]
                    assign("current.best", current.best, envir = envmh)
                    assign("current.best.GP", GP, envir = envmh)
                    X <- FUN.I(envmh = envmh, X = X)
                }
            }
        }
    }
    if (!is.null(FUN.EXIT)) {
        FUN.EXIT(envmh = envmh, X = X)
    }
    assign("GP", GP, pos = globalenv())
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
