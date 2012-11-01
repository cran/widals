

Z <- read.table("http://www.stat.ucla.edu/~davezes/code_and_data/O3.tsv", sep = "\t", stringsAsFactors = FALSE)

locs <- read.table("http://www.stat.ucla.edu/~davezes/code_and_data/locs.tsv", sep = "\t", stringsAsFactors = FALSE)

if(FALSE) {
load("http://www.stat.ucla.edu/~davezes/code_and_data/O3stoFull_129.robj")

load("~/Desktop/O3stoFull_129.robj")

names(O3.sto.f129)

dim(O3.sto.f129$Z)


O3.sto.f129$date

O3.sto.f129$na.mask
}

save(Z, locs, file="~/Desktop/O3.RData")


library(sp)
library(maps)
data(O3)


min(locs[ ,2])
max(locs[ ,2])

min(locs[ ,1])
max(locs[ ,1])

xmp <- map("state", "california", plot=FALSE)
plot(xmp, type="n")
lines(xmp)

points(locs)


xlat.dom <- seq(min(locs[ ,2])-0.5, max(locs[ ,2])+2, length=84)
xlon.dom <- seq(min(locs[ ,1])-0.9, max(locs[ ,1])+0.5, length=70)
locs0 <- cbind(rep(xlat.dom, length(xlon.dom)), rep(xlon.dom, each=length(xlat.dom)))[ , c(2,1)]
points(locs0)

#xint.mask <- point.in.polygon(locs0[ ,1], locs0[ ,2], xmp$x, xmp$y)
#xmask <- rep(FALSE, nrow(locs0))
#xmask[xint.mask==2 ] <- TRUE


smask <- map.where(database = "state", locs0[ ,1], locs0[ ,2])
xmask <- rep(FALSE, nrow(locs0))
xmask[smask=="california" ] <- TRUE ; xmask[ is.na(xmask) ] <- TRUE
sum(xmask)
locs0x <- locs0[xmask, ]

xmp <- map("state", "california", plot=FALSE)
plot(xmp, type="n")
lines(xmp)
points(locs0x)


#zzlocs <- locs #### edit point
zzlocs <- locs0x

xelevs <- NULL
for(i in 1:nrow(zzlocs)) {
    
    this.lat <- zzlocs[i, 2]
    this.lon <- zzlocs[i, 1]
    xsource <- paste( "http://maps.googleapis.com/maps/api/elevation/json?locations=", this.lat, ",", this.lon, "&sensor=false",    sep="" )
    aa <- read.table( xsource, sep="\n", stringsAsFactors=FALSE )
    bb <- strsplit( aa[ 4, ] , ":", fixed=TRUE )
    gg <- bb[[1]][2] ; gg
    gg <- gsub(" ", "", gg, fixed=TRUE)
    gg <- gsub(",", "", gg, fixed=TRUE)
    xelevs[i] <- gg
    cat(i, gg, "\n")
    
}

#write.table( xelevs, "~/Desktop/Testelevs.csv", row.names=FALSE, quote=FALSE, sep=",")

xelevs <- round(as.numeric(xelevs)*10) / 10
CAelevs <- xelevs


save(CAelevs, locs0x, file="~/Desktop/CAelevs.robj")



save(Z, locs, xelevs, file="~/Files/Creations/R/widals/Data/O3.RData")








xelevs <- rep(NA, nrow(zzlocs))
for(i in 2320:nrow(zzlocs)) {
    
    this.lat <- zzlocs[i, 2]
    this.lon <- zzlocs[i, 1]
    xsource <- paste( "http://maps.googleapis.com/maps/api/elevation/json?locations=", this.lat, ",", this.lon, "&sensor=false",    sep="" )
    aa <- read.table( xsource, sep="\n", stringsAsFactors=FALSE )
    bb <- strsplit( aa[ 4, ] , ":", fixed=TRUE )
    gg <- bb[[1]][2] ; gg
    gg <- gsub(" ", "", gg, fixed=TRUE)
    gg <- gsub(",", "", gg, fixed=TRUE)
    xelevs[i] <- gg
    cat(i, gg, "\n")
    
}

xxxxx <- c(CAelevs, round(as.numeric(xelevs[ 2322:length(xelevs) ])*10) / 10 )

CAelevs <- xxxxx

library(widals)
library(maps)
data(O3)

load(file="~/Desktop/CAelevs.robj")

save(Z, locs, xelevs, CAelevs, locs0x, file="~/Files/Creations/R/widals/Data/O3.RData")


xmp <- map("state", "california", plot=FALSE)
plot(xmp, type="n")
lines(xmp)
points(locs0x[ 1:length(CAelevs), ], cex=CAelevs/1000)





locs <- read.table("http://www.stat.ucla.edu/~davezes/code_and_data/locs.tsv", sep = "\t", stringsAsFactors = FALSE)








options(stringsAsFactors=FALSE)

k.cpus <- 5 #### set the number of cpus for snowfall


library(widals)
tau <- 210
n.all <- 300

set.seed(77777)
locs.all <- cbind(runif(n.all), runif(n.all))

D.mx <- distance(locs.all, locs.all, FALSE)

Q <- 0.03*exp(-2*D.mx)
F <- 0.99
R <- diag(1, n.all)
beta0 <- rep(0, n.all)
H <- diag(1, n.all)

xsssim <- SS.sim(F, H, Q, R, length.out=tau, beta0=beta0)
Y1 <- xsssim$Y



Hst.ss <- list()
for(tt in 1:tau) {
Hst.ss[[tt]] <- cbind( rep(sin(tt*2*pi/tau), n.all), rep(cos(tt*2*pi/tau), n.all) )
colnames(Hst.ss[[tt]]) <- c("sinet", "cosinet")
}
Ht.original <- cbind( sin((1:tau)*2*pi/tau), cos((1:tau)*2*pi/tau) )

Q2 <- diag(0.03, ncol(Hst.ss[[1]]))
F2 <- 0.99
beta20 <- rep(0, ncol(Hst.ss[[1]]))
R2 <- 1*exp(-3*D.mx) + diag(0.001, n.all)

xsssim2 <- SS.sim.tv(F2, Hst.ss, Q2, R, length.out=tau, beta0=beta20)
Z2 <- xsssim2$Z


Z.all <- Y1 + Z2



z.min <- min(Z.all)
for(tt in 1:tau) {
plot(locs.all, cex=(Z.all[ tt, ]-z.min)*0.3,   main=tt)
Sys.sleep(0.1)
}





train.rng <- 30:tau
test.rng <- train.rng


xgeodesic <- FALSE
Z <- Z.all
locs <- locs.all
n <- n.all
Ht <- Ht.original
Hs <- matrix(1, nrow=n)
Hst.ls <- NULL


rm.ndx <- create.rm.ndx.ls(n, 14)
k.glob <- 10
run.parallel <- TRUE

FUN.source <- fun.load.widals.a

d.alpha.lower.limit <- 0
rho.upper.limit <- 100
rgr.lower.limit <- 10^(-7)

GP <- c(1/10, 1, 0.01, 3, 1)
cv <- 2
lags <- c(0)
b.lag <- 0

sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))
ltco <- -10
stnd.d <- TRUE

sfInit(TRUE, k.cpus)
FUN.source()

set.seed(99999)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx,
k.glob, k.loc.coef=7, X = NULL)
sfStop()

### 1.102488






library(LatticeKrig)

xxb <- 7
yyb <- 7
center <- as.matrix(expand.grid( seq( 0, 1, length=xxb), seq( 0, 1, length=yyb)))

###### run together
ffunit <- function(x) { return( (x-min(x)) / (max(x)-min(x)) ) }
locs.unit <- apply(locs, 2, ffunit)
locs.unit <- matrix(as.vector(locs.unit), ncol=2)
###### run togther

xPHI <- Radial.basis(as.matrix(locs.unit), center, 0.5)
Hs.lkrig <- matrix(NA, nrow(locs), xxb*yyb)
for(i in 1:nrow(locs)) {
    Hs.lkrig[i, ] <- xPHI[i]
}




XA.Hs <- Hs
XA.Ht <- Ht.original
XA.Hst.ls <- NULL

Hst.sumup(XA.Hst.ls, XA.Hs, XA.Ht)


XB.Hs <- 10*Hs.lkrig
XB.Ht <- NULL
XB.Hst.ls <- NULL

Hst.sumup(XB.Hst.ls, XB.Hs, XB.Ht)



GP <- c(1/10, 1, 1/10, 1,      5, 3, 1)
sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))

rm.ndx <- create.rm.ndx.ls(n, 14)
cv <- 2
lags <- c(0)
b.lag <- 0

d.alpha.lower.limit <- 0
rho.upper.limit <- 100
rgr.lower.limit <- 10^(-7)
rgr.upper.limit <- 100

source("/Users/dzes/Desktop/newfun.R")

sfInit(TRUE, k.cpus)
FUN.source <- fun.load.widals.ab
FUN.source()

set.seed(9999)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx,
k.glob, k.loc.coef=7, X = NULL)
sfStop()

#### 












