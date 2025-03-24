widals.predict <-
function (Z, Hs, Ht, Hst.ls, locs, lags, b.lag, Hs0, Hst0.ls, 
    locs0, geodesic = FALSE, wrap.around = NULL, GP, stnd.d = FALSE, 
    ltco = -16) 
{
    tau <- nrow(Z)
    n <- nrow(locs)
    k <- length(lags)
    n0 <- nrow(locs0)
    rho <- GP[1]
    reg <- GP[2]
    alpha <- GP[3]
    beta <- GP[4]
    flatten <- GP[5]
    locs0.3D <- cbind(locs0, rep(0, n0))
    locs.long.3D <- cbind(rep(locs[, 1], k), rep(locs[, 2], k), 
        beta * rep(lags, each = n))
    z.lags.vec <- rep(lags, each = n)
    ALS <- H.als.b(Z = Z, Hs = Hs, Ht = Ht, Hst.ls = Hst.ls, 
        rho = rho, reg = reg, b.lag = b.lag, Hs0 = Hs0, Ht0 = Ht, 
        Hst0.ls = Hst0.ls)
    Y.als <- ALS$Z.hat
    dim(Y.als)
    Y0.als <- ALS$Z0.hat
    dim(Y0.als)
    rm(ALS)
    Z.delta <- Z - Y.als
    Z.delta <- Z.clean.up(Z.delta)
    Z.adj <- crispify(locs1 = locs0.3D, locs2 = locs.long.3D, 
        Z.delta = Z.delta, z.lags.vec = z.lags.vec, geodesic = geodesic, 
        alpha = alpha, flatten = flatten, self.refs = c(-1), 
        lags = lags, stnd.d = stnd.d, log10cutoff = ltco)
    Z0.wid <- Y0.als + Z.adj
    return(Z0.wid)
}
