fun.load.hals.a <-
function() {
    
    
    
    if( run.parallel ) {
        
        sfExport("Z", "Hs", "Ht", "Hst.ls", "b.lag", "train.rng", "test.rng")
        
        suppressWarnings(sfLibrary(widals))
        
    }
    
    p.ndx.ls <- list( c(1,2) )
    assign( "p.ndx.ls", p.ndx.ls, pos=globalenv() )
    
    f.d <- list( dlog.norm, dlog.norm, dlog.norm, dlog.norm, dlog.norm )
    assign( "f.d", f.d, pos=globalenv() )
    
    
    FUN.MH <- function(jj, GP.mx, X) {
        
        Z.als <- Hals.snow( jj, Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, b.lag=b.lag, GP.mx=GP.mx )
        
        resids <- Z - Z.als
        our.cost <- sqrt( mean( resids[ train.rng, ]^2 ) )
        
        return( our.cost )
    }

    assign( "FUN.MH", FUN.MH, pos=globalenv() )
    
    #FUN.GP <- NULL
    
    
    FUN.GP <- function(GP.mx) {
        GP.mx[ GP.mx[ , 1] > rho.upper.limit, 1 ] <- rho.upper.limit
        GP.mx[ GP.mx[ , 2] < rgr.lower.limit, 2 ] <- rgr.lower.limit
        return(GP.mx)
        
    }
    assign( "FUN.GP", FUN.GP, pos=globalenv() )
    
    
    FUN.I <- function(envmh, X) {
        cat( "Improvement ---> ", envmh$current.best, " ---- " , envmh$GP, "\n" )
        
    }
    assign( "FUN.I", FUN.I, pos=globalenv() )
    
    
    
    FUN.EXIT <- function(envmh, X) {
        
        GP.mx <- matrix(envmh$GP, 1, length(envmh$GP))
        
        Z.als <- Hals.snow( 1, Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, b.lag=b.lag, GP.mx=GP.mx )
        
        resids <- Z - Z.als
        our.cost <- sqrt( mean( resids[ test.rng, ]^2 ) )
        
        cat( envmh$GP, " -- ", our.cost,   "\n" )
        
        assign( "Z.als", Z.als, pos=globalenv() )
        assign( "our.cost", our.cost, pos=globalenv() )
        assign( "GP", envmh$GP, pos=globalenv() )
        cat( paste( "GP <- c(", paste(format(GP,digits=5), collapse=", "), ") ### ", format(our.cost, width=6), "\n", sep="" ) )
        
    }
    
    assign( "FUN.EXIT", FUN.EXIT, pos=globalenv() )
    
}
