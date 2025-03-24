\name{fun.load}
\alias{fun.load}
\alias{fun.load.hals.a}
\alias{fun.load.hals.fill}
\alias{fun.load.widals.a}
\alias{fun.load.widals.fill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic Search Helper Functions
}
\description{
Functions that assign values and functions needed by \code{\link{MSS.snow}}
}
\usage{
fun.load.hals.a()
fun.load.hals.fill()
fun.load.widals.a()
fun.load.widals.fill()
}
%- maybe also 'usage' for other objects documented here.
\details{
Please see \code{\link{MSS.snow}} and examples.
}
\value{
Nothing.  The central role of these functions is the creation of four functions required by \code{\link{MSS.snow}}: \code{FUN.MH}, \code{FUN.GP}, \code{FUN.I}, and \code{FUN.EXIT}.  These four functions are assigned to the Global Environment.  This \code{\link{fun.load}} suite of functions also passes needed objects (out-of-scope) to \code{\link[snowfall]{snowfall}} threads if the global user-made variable \code{run.parallel} is set to \code{TRUE}.
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
\code{\link{MSS.snow}}
}


\examples{



### Here's an itty bitty example:
### we use stochastic search to find the minimum number in a vector
### GP isn't used here, and hence neither are p.ndx.ls nor f.d
### however, we still need to create them since MSS.snow requires their existence

\dontrun{

fun.load.simpleExample <- function() {

   if( run.parallel ) {
         sfExport("xx")
    }
    
    p.ndx.ls <- list( c(1) )
    
    p.ndx.ls <<- p.ndx.ls
    
    f.d <- list( dlog.norm )

    f.d <<- f.d
    
    FUN.MH <- function(jj, GP.mx, X) {
        our.cost <- sample(xx, 1)
    }

    FUN.MH <<- FUN.MH
    
    
    FUN.GP <- NULL
    FUN.GP <<- FUN.GP
    
    
    FUN.I <- function(envmh, X) {
        cat( "Hello, I have found an even smaller number in xx ---> ", envmh$current.best, "\n" )
    }
    FUN.I <<- FUN.I
    
    FUN.EXIT <- function(envmh, X) {
        cat( "Done",   "\n" )
    }

    FUN.EXIT <<- FUN.EXIT
    
}

xx <- 1:600

GP <- c(1)

run.parallel <- TRUE
sfInit(TRUE, 2)

MH.source <- fun.load.simpleExample
MH.source()

MSS.snow(MH.source, Inf, p.ndx.ls, f.d, matrix(1, nrow=28), 28, 7)
sfStop()




### Here's another itty bitty example:
### we use stochastic search to find the mean of a vector
### i.e., the argmin? of sum ( x - ? )^2

fun.load.simpleExample2 <- function() {

   if( run.parallel ) {
         sfExport("xx")
    }
    
    p.ndx.ls <- list( c(1) )
    p.ndx.ls <<- p.ndx.ls
    
    f.d <- list( unif.mh )
    f.d <<- f.d
    
    FUN.MH <- function(jj, GP.mx, X) {
        our.cost <- sum( ( xx - GP.mx[jj, 1] )^2 )
        return(our.cost)
    }
    FUN.MH <<- FUN.MH
    
    FUN.GP <- NULL
    FUN.GP <<-  FUN.GP
    
    FUN.I <- function(envmh, X) {
        cat( "Improvement ---> ", envmh$current.best, " ---- " , envmh$GP, "\n" )
    }
    FUN.I <<- FUN.I
    
    FUN.EXIT <- function(envmh, X) {
        our.cost <- envmh$current.best
        GP <- envmh$GP
        cat( "Done",   "\n" )
        cat( envmh$GP, our.cost, "\n" )
    }
    FUN.EXIT <<- FUN.EXIT
    
}

##set.seed(99999)
xx <- rnorm(300, 5, 10)

GP <- c(1)

run.parallel <- TRUE
sfInit(TRUE, 2)

MH.source <- fun.load.simpleExample2
MH.source()

MSS.snow(MH.source, Inf, p.ndx.ls, f.d, matrix(1/10, nrow=140, ncol=length(GP)), 140, 14)
sfStop()

##### in fact:
mean(xx)

}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


