fun.load.hals.a <-
function () 
{
    run.parallel <- run.parallel
    if (run.parallel) {
        sfExport("Z", "Hs", "Ht", "Hst.ls", "b.lag", "train.rng", 
            "test.rng")
        suppressWarnings(sfLibrary("widals", character.only = TRUE))
    }
    p.ndx.ls <- list(c(1, 2))
    p.ndx.ls <<- p.ndx.ls
    f.d <- list(dlog.norm, dlog.norm, dlog.norm, dlog.norm, dlog.norm)
    f.d <<- f.d
    FUN.MH <- function(jj, GP.mx, X) {
        Z <- Z
        Hs <- Hs
        Ht <- Ht
        Hst.ls <- Hst.ls
        b.lag <- b.lag
        train.rng <- train.rng
        Z.als <- Hals.snow(jj, Z = Z, Hs = Hs, Ht = Ht, Hst.ls = Hst.ls, 
            b.lag = b.lag, GP.mx = GP.mx)
        resids <- Z - Z.als
        our.cost <- sqrt(mean(resids[train.rng, ]^2))
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
        cat("Improvement ---> ", envmh$current.best, " ---- ", 
            envmh$GP, "\n")
    }
    FUN.I <<- FUN.I
    FUN.EXIT <- function(envmh, X) {
        Z <- Z
        Hs <- Hs
        Ht <- Ht
        Hst.ls <- Hst.ls
        b.lag <- b.lag
        train.rng <- train.rng
        test.rng <- test.rng
        GP.mx <- matrix(envmh$GP, 1, length(envmh$GP))
        Z.als <- Hals.snow(1, Z = Z, Hs = Hs, Ht = Ht, Hst.ls = Hst.ls, 
            b.lag = b.lag, GP.mx = GP.mx)
        resids <- Z - Z.als
        our.cost <- sqrt(mean(resids[test.rng, ]^2))
        cat(envmh$GP, " -- ", our.cost, "\n")
        Z.als <<- Z.als
        our.cost <<- our.cost
        GP <- envmh$GP
        GP <<- GP
        cat(paste("GP <- c(", paste(format(GP, digits = 5), collapse = ", "), 
            ") ### ", format(our.cost, width = 6), "\n", sep = ""))
    }
    FUN.EXIT <<- FUN.EXIT
}
