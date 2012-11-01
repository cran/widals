y.unalign <-
function(Y.sub, lags) {
	tau0 <- nrow(Y.sub)
	low.ndx <- min( min(lags), 0 )
	top.ndx <- max( 0, max(lags) )
	tau <- tau0 + (top.ndx - low.ndx)
	Y <- matrix( NA, tau, ncol(Y.sub) )
	Y[ (1-low.ndx):(tau0-low.ndx), ] <- Y.sub
	return(Y)
}
