fun.load.widals.a <-
function() {
    
    
    
    if( run.parallel ) {
        
        sfExport("Z", "Hs", "Ht", "Hst.ls", "locs", "lags", "b.lag", "cv", "rm.ndx", "train.rng", "test.rng", "xgeodesic", "ltco", "stnd.d")
        
        suppressWarnings(sfLibrary(widals))
        
    }
    
    if( length(lags) == 1 & lags[1] == 0 ) {
        p.ndx.ls <- list( c(1,2), c(3,5) )
    } else {
        p.ndx.ls <- list( c(1,2), c(3,4,5) )
    }
    assign( "p.ndx.ls", p.ndx.ls, pos=globalenv() )
    
    f.d <- list( dlog.norm, dlog.norm, dlog.norm, dlog.norm, dlog.norm )
    assign( "f.d", f.d, pos=globalenv() )
    
    
    FUN.MH <- function(jj, GP.mx, X) {
        
        Z.wid <- widals.snow(jj, rm.ndx=rm.ndx, Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, locs=locs, lags=lags, b.lag=b.lag, cv=cv, geodesic=xgeodesic,
        wrap.around=NULL, GP.mx, stnd.d=stnd.d, ltco=ltco)
        
        if( min(Z, na.rm=TRUE) >= 0 ) { Z.wid[ Z.wid < 0 ] <- 0 } ############ DZ EDIT
        
        Z.wid <- Z.clean.up(Z.wid)
        
        resids <- Z[ , unlist(rm.ndx)] - Z.wid[ , unlist(rm.ndx)]
        our.cost <- sqrt( mean( resids[ train.rng, ]^2 ) )
        
        if( is.nan(our.cost) ) { our.cost <- Inf }
        
        return( our.cost )
    }
    assign( "FUN.MH", FUN.MH, pos=globalenv() )
    
    #FUN.GP <- NULL
    
    
    FUN.GP <- function(GP.mx) {
        GP.mx[ GP.mx[ , 1] > rho.upper.limit, 1 ] <- rho.upper.limit
        GP.mx[ GP.mx[ , 2] < rgr.lower.limit, 2 ] <- rgr.lower.limit
        GP.mx[ GP.mx[ , 3] < d.alpha.lower.limit, 3 ] <- d.alpha.lower.limit
        xperm <- order(GP.mx[ , 3,  drop=FALSE])
        GP.mx <- GP.mx[ xperm,  ,  drop=FALSE]
        
        return(GP.mx)
        
    }
    assign( "FUN.GP", FUN.GP, pos=globalenv() )
    
    
    FUN.I <- function(envmh, X) {
        cat( "Improvement ---> ", envmh$current.best, " ---- " , envmh$GP, "\n" )
        
    }
    assign( "FUN.I", FUN.I, pos=globalenv() )
    
    
    
    FUN.EXIT <- function(envmh, X) {
        
        GP.mx <- matrix(envmh$GP, 1, length(envmh$GP))
        Z.wid <- widals.snow(1, rm.ndx=rm.ndx, Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, locs=locs, lags=lags, b.lag=b.lag, cv=cv, geodesic=xgeodesic,
        wrap.around=NULL, GP.mx, stnd.d=stnd.d, ltco=ltco)
        
        if( min(Z, na.rm=TRUE) >= 0 ) { Z.wid[ Z.wid < 0 ] <- 0 } ############ DZ EDIT
        
        assign( "Z.wid", Z.wid, envir=globalenv() )
        
        Z.wid <- Z.clean.up(Z.wid)
        
        resids <- Z[ , unlist(rm.ndx)] - Z.wid[ , unlist(rm.ndx)]
        our.cost <- sqrt( mean( resids[ test.rng, ]^2 ) )
        
        if( is.nan(our.cost) ) { our.cost <- Inf }
        
        cat( envmh$GP, " -- ", our.cost,   "\n" )
        
        assign( "our.cost", our.cost, pos=globalenv() )
        assign( "GP", envmh$GP, pos=globalenv() )
        cat( paste( "GP <- c(", paste(format(GP,digits=5), collapse=", "), ") ### ", format(our.cost, width=6), "\n", sep="" ) )
        
    }
    
    assign( "FUN.EXIT", FUN.EXIT, pos=globalenv() )
    
}
