fun.load.hals.fill <-
function () 
{
    run.parallel <- run.parallel
    if (run.parallel) {
        sfExport("Z", "Hs", "Ht", "Hst.ls", "b.lag", "train.rng", 
            "test.rng", "Z.na")
        suppressWarnings(sfLibrary("widals", character.only = TRUE))
    }
    p.ndx.ls <- list(c(1, 2))
    p.ndx.ls <<- p.ndx.ls
    f.d <- list(dlog.norm, dlog.norm, dlog.norm, dlog.norm, dlog.norm)
    f.d <<- f.d
    FUN.MH <- function(jj, GP.mx, X) {
        Z.na <- Z.na
        Z <- Z
        Hs <- Hs
        Ht <- Ht
        Hst.ls <- Hst.ls
        b.lag <- b.lag
        train.rng <- train.rng
        Z.als <- Hals.snow(jj, Z = X$Z.fill, Hs = Hs, Ht = Ht, 
            Hst.ls = Hst.ls, b.lag = b.lag, GP.mx = GP.mx)
        resids <- (Z - Z.als)[train.rng, ]
        our.cost <- sqrt(mean(resids[!Z.na[train.rng, ]]^2))
        return(our.cost)
    }
    FUN.MH <<- FUN.MH
    FUN.GP <- function(GP.mx) {
        rho.upper.limit <- rho.upper.limit
        rgr.lower.limit <- rgr.lower.limit
        GP.mx[GP.mx[, 1] > rho.upper.limit, 1] <- rho.upper.limit
        GP.mx[GP.mx[, 2] < rgr.lower.limit, 2] <- rgr.lower.limit
        return(GP.mx)
    }
    FUN.GP <<- FUN.GP
    FUN.I <- function(envmh, X) {
        Hs <- Hs
        Ht <- Ht
        Hst.ls <- Hst.ls
        b.lag <- b.lag
        Z.na <- Z.na
        GP.mx <- matrix(envmh$GP, 1, length(envmh$GP))
        Z.als <- Hals.snow(1, Z = X$Z.fill, Hs = Hs, Ht = Ht, 
            Hst.ls = Hst.ls, b.lag = b.lag, GP.mx = GP.mx)
        Z.als <<- Z.als
        Z.als <- Z.clean.up(Z.als)
        X$Z.fill[Z.na] <- Z.als[Z.na]
        X$Z.fill[is.na(X$Z.fill)] <- mean(X$Z.fill, na.rm = TRUE)
        cat("Improvement ---> ", envmh$current.best, " ---- ", 
            envmh$GP, "\n")
        return(X)
    }
    FUN.I <<- FUN.I
    FUN.EXIT <- function(envmh, X) {
        our.cost <- envmh$current.best
        our.cost <<- our.cost
        Z.fill <- X$Z.fill
        Z.fill <<- Z.fill
        GP <- envmh$GP
        GP <<- GP
        cat(paste("GP <- c(", paste(format(envmh$GP, digits = 5), 
            collapse = ", "), ") ### ", format(our.cost, width = 6), 
            "\n", sep = ""))
    }
    FUN.EXIT <<- FUN.EXIT
}
