fun.load.widals.fill <-
function() {
    
    lags <- lags
    run.parallel <- run.parallel

    
    if( run.parallel ) {
        
        sfExport("Z", "Hs", "Ht", "Hst.ls", "locs", "lags", "b.lag", "cv", "rm.ndx", "train.rng", "test.rng", "xgeodesic", "ltco", "stnd.d", "Z.na")
        
        suppressWarnings(sfLibrary("widals", character.only=TRUE))
        
    }
    
    
    if( length(lags) == 1 & lags[1] == 0 ) {
        p.ndx.ls <- list( c(1,2), c(3,5) )
    } else {
        p.ndx.ls <- list( c(1,2), c(3,4,5) )
    }
    ## assign( "p.ndx.ls", p.ndx.ls, pos=globalenv() )
    p.ndx.ls <<- p.ndx.ls
    
    f.d <- list( dlog.norm, dlog.norm, dlog.norm, dlog.norm, dlog.norm )
    ## assign( "f.d", f.d, pos=globalenv() )
    f.d <<- f.d
    
    FUN.MH <- function(jj, GP.mx, X) {
        
        rm.ndx <- rm.ndx
        Hs <- Hs
        Ht <- Ht
        Hst.ls <- Hst.ls
        locs <- locs
        lags <- lags
        b.lag <- b.lag
        cv <- cv
        xgeodesic <- xgeodesic
        stnd.d <- stnd.d
        ltco <- ltco
        Z <- Z
        train.rng <- train.rng
        Z.na <- Z.na
        
        Z.wid <- widals.snow(jj, rm.ndx=rm.ndx, Z=X$Z.fill, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, locs=locs, lags=lags, b.lag=b.lag,
        cv=cv, geodesic=xgeodesic, wrap.around=NULL, GP.mx=GP.mx, stnd.d=stnd.d, ltco=ltco )
        
        
        ##if( min(Z, na.rm=TRUE) >= 0 ) { Z.wid[ Z.wid < 0 ] <- 0 } ############ DZ EDIT
        
        resids <- (Z - Z.wid)[ train.rng, unlist(rm.ndx) ]
        our.cost <- sqrt( mean( resids[ !Z.na[ train.rng, unlist(rm.ndx) ] ]^2 ) )
        
        return( our.cost )
    }
    ## assign( "FUN.MH", FUN.MH, pos=globalenv() )
    FUN.MH <<- FUN.MH
    
    #FUN.GP <- NULL
    
    
    FUN.GP <- function(GP.mx) {
        
        rho.upper.limit <- rho.upper.limit
        rgr.lower.limit <- rgr.lower.limit
        d.alpha.lower.limit <- d.alpha.lower.limit
        
        
        GP.mx[ GP.mx[ , 1] > rho.upper.limit, 1 ] <- rho.upper.limit
        GP.mx[ GP.mx[ , 2] < rgr.lower.limit, 2 ] <- rgr.lower.limit
        GP.mx[ GP.mx[ , 3] < d.alpha.lower.limit, 3 ] <- d.alpha.lower.limit
        xperm <- order(GP.mx[ , 3,  drop=FALSE])
        GP.mx <- GP.mx[ xperm,  ,  drop=FALSE]
        
        return(GP.mx)
        
    }
    ## assign( "FUN.GP", FUN.GP, pos=globalenv() )
    FUN.GP <<- FUN.GP
    
    FUN.I <- function(envmh, X) {
        
        rm.ndx <- rm.ndx
        
        Hs <- Hs
        Ht <- Ht
        Hst.ls <- Hst.ls
        locs <- locs
        lags <- lags
        
        b.lag <- b.lag
        cv <- cv
        xgeodesic <- xgeodesic
        stnd.d <- stnd.d
        ltco <- ltco
        Z.na <- Z.na
        

        
        GP.mx <- matrix(envmh$GP, 1, length(envmh$GP))
        Z.wid <- widals.snow(1, rm.ndx=rm.ndx, Z=X$Z.fill, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, locs=locs, lags=lags, b.lag=b.lag,
        cv=cv, geodesic=xgeodesic, wrap.around=NULL, GP.mx=GP.mx, stnd.d=stnd.d, ltco=ltco )
        
        ### if( min(Z, na.rm=TRUE) >= 0 ) { Z.wid[ Z.wid < 0 ] <- 0 } ############ DZ EDIT
        
        ## assign( "Z.wid", Z.wid, envir=globalenv() )
        Z.wid <<- Z.wid
        Z.wid <- Z.clean.up(Z.wid)
        
        
        X$Z.fill[ Z.na ]  <- Z.wid[ Z.na ]
        X$Z.fill[ is.na(X$Z.fill) ] <- mean(  X$Z.fill, na.rm=TRUE )
        
        
        cat( "Improvement ---> ", envmh$current.best, " ---- " , envmh$GP, "\n" )
        return(X)
    }
    ## assign( "FUN.I", FUN.I, pos=globalenv() )
    FUN.I <<- FUN.I
    
    FUN.EXIT <- function(envmh, X) {
        
        our.cost <- envmh$current.best
        ## assign( "our.cost", our.cost, pos=globalenv() )
        ## assign( "Z.fill", X$Z.fill, envir=globalenv() )
        ## assign( "GP", envmh$GP, pos=globalenv() )
        
        our.cost <<- our.cost
        Z.fill <- X$Z.fill
        Z.fill <<- Z.fill
        GP <- envmh$GP
        GP <<- GP
        
        cat( paste( "GP <- c(", paste(format(envmh$GP,digits=5), collapse=", "), ") ### ", format(our.cost, width=6), "\n", sep="" ) )
        
        
    }
    ## assign( "FUN.EXIT", FUN.EXIT, pos=globalenv() )
    FUN.EXIT <<- FUN.EXIT
    
}
