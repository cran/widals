als.prepare <-
function(X, Z, lags, Z.na, tt.rng, rgr=NULL, assignToEnvironment=NULL) {
	yqq <- (tt.rng[1] + min(lags, 0))
	y.tt.rng <- I( yqq:(yqq+length(tt.rng)-1) )
    
	if( is.null(X) ) {
		xalgn <- ar.align( Z, lags, Z.na )
	} else {
		xalgn <- xy.align( X, Z, lags, Z.na )
	}
	X.sub <- cbind( xalgn$X.sub, rgr[ xalgn$y.sub.ndx, , drop=FALSE ] )
	if( is.null(assignToEnvironment ) ) {
		return( list(X.sub=X.sub, Y.sub=xalgn$Y.sub, Y.na=xalgn$Y.na, y.tt.rng=y.tt.rng, y.sub.ndx=xalgn$y.sub.ndx) )
	} else {
		assign( "X.sub", X.sub, pos=assignToEnvironment )
		assign( "Y.sub", xalgn$Y.sub, pos=assignToEnvironment )
		assign( "Y.na", xalgn$Y.na, pos=assignToEnvironment )
		assign( "y.tt.rng", y.tt.rng, pos=assignToEnvironment )
		assign( "y.sub.ndx", xalgn$y.sub.ndx, pos=assignToEnvironment )
	}
    
}
