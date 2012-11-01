




options(stringsAsFactors=FALSE)

k.cpus <- 5 #### set the number of cpus for snowfall

library(widals)
data(O3)

Z.all <- as.matrix(O3$Z)[366:730, ]
locs.all <- O3$locs[ , c(2,1)]
hsa.all <- O3$helevs/500

xdate <- rownames(Z.all)

tau <- nrow(Z.all)
n.all <- ncol(Z.all)

xgeodesic <- TRUE


Z <- Z.all
locs <- locs.all
n <- n.all

dateDate <- strptime(xdate, "%Y%m%d")
doy <- as.integer(format(dateDate, "%j"))

Ht <- cbind( sin(2*pi*doy/365), cos(2*pi*doy/365) )

Hs.all <- cbind(matrix(1, nrow=n.all), hsa.all)

Hisa.ls <- H.Earth.solar(locs[ , 2], locs[ , 1], dateDate)

Hst.ls.all2 <- list()
for(tt in 1:tau) {
    Hst.ls.all2[[tt]] <- cbind(Hisa.ls[[tt]], Hisa.ls[[tt]]*hsa.all)
    colnames(Hst.ls.all2[[tt]]) <- c("ISA", "ISAxElev")
}
Hst.ls <- Hst.ls.all2

Hs <- Hs.all

Ht.original <- Ht


train.rng <- 30:tau
test.rng <- train.rng

rm.ndx <- 1:n
k.glob <- 10
run.parallel <- TRUE

FUN.source <- fun.load.widals.a

d.alpha.lower.limit <- 0
rho.upper.limit <- 100
rgr.lower.limit <- 10^(-7)

GP <- c(1/10, 1, 0.01, 3, 1)
cv <- -2
lags <- c(0)
b.lag <- -1

sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))
ltco <- -10
stnd.d <- TRUE

sfInit(TRUE, k.cpus)
FUN.source()

set.seed(99999)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx,
k.glob, k.loc.coef=7, X = NULL)
sfStop()

#### 11.90536




rm.ndx <- create.rm.ndx.ls(n, 14)
k.glob <- 10

FUN.source <- fun.load.widals.a

GP <- c(1/10, 1, 0.01, 3, 1)
cv <- 2
lags <- c(0)
b.lag <- 0

sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))
ltco <- -10

sfInit(TRUE, k.cpus)
FUN.source()

set.seed(99999)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx,
k.glob, k.loc.coef=7, X = NULL)
sfStop()

#### 12.11686




library(LatticeKrig)

xxb <- 4
yyb <- 4
center <- as.matrix(expand.grid( seq( 0, 1, length=xxb), seq( 0, 1, length=yyb)))

locs.unit <- as.matrix(t( (t(locs) - apply(locs, 2, min)) / apply(locs, 2, function(x) { return(max(x)-min(x)) } ) ))
locs.unit <- matrix(as.vector(locs.unit), ncol=2)

xPHI <- Radial.basis(as.matrix(locs.unit), center, 0.5)
Hs.lkrig <- matrix(NA, nrow(locs), xxb*yyb)
for(i in 1:nrow(locs)) {
Hs.lkrig[i, ] <- xPHI[i]
}
Hs.lkrig <- Hs.lkrig[ , -which( apply(Hs.lkrig, 2, sum) < 0.1 ), drop=FALSE ]








Hs <- 10*Hs.lkrig

Hs <- cbind(rep(1,n), hsa.all, 10*Hs.lkrig ) ### 12.1124

Hs <- cbind(rep(1,n), hsa.all)

diag(Hst.sumup(Hst.ls, Hs, Ht))


GP <- c(1/10, 1, 5, 3, 1)
rm.ndx <- create.rm.ndx.ls(n, 14)
cv <- 2
lags <- c(0)
b.lag <- 0

sfInit(TRUE, k.cpus)
FUN.source()

MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx,
k.glob, k.loc.coef=7, X = NULL)
sfStop()








diag(Hst.sumup(Hst.ls, Hs, Ht))









XA.Hs <- cbind(rep(1,n), hsa.all)
XA.Ht <- Ht.original
XA.Hst.ls <- Hst.ls

Hst.sumup(XA.Hst.ls, XA.Hs, XA.Ht)


XB.Hs <- 10*Hs.lkrig
XB.Ht <- NULL
XB.Hst.ls <- NULL

Hst.sumup(XB.Hst.ls, XB.Hs, XB.Ht)



GP <- c(1/10, 1, 1/10, 1,      5, 3, 1)
sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))

rm.ndx <- I(1:n)
cv <- -2
lags <- c(0)
b.lag <- -1

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

#### 11.5234




GP <- c(1/10, 1, 1/10, 1,      5, 3, 1)
sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))

rm.ndx <- create.rm.ndx.ls(n, 14)
cv <- 2
lags <- c(0)
b.lag <- 0

d.alpha.lower.limit <- 0
rho.upper.limit <- 100
rgr.lower.limit <- 10^(-7)
rgr.upper.limit <- 300

source("/Users/dzes/Desktop/newfun.R")
sfInit(TRUE, k.cpus)
FUN.source <- fun.load.widals.ab
FUN.source()

set.seed(9999)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx,
k.glob, k.loc.coef=9, X = NULL)
sfStop()






##################################

matrix(Radial.basis(as.matrix(locs.unit), center, 0.5), nrow=nrow(locs), byrow=TRUE)






center <- as.matrix(expand.grid( seq( min(locs[,1]), max(locs[,1]), length=3), seq( min(locs[,2]), max(locs[,2]), length=3)))

locs.unit <- as.matrix(t( (t(locs) - apply(locs, 2, min)) / apply(locs, 2, function(x) { return(max(x)-min(x)) } ) ))
locs.unit <- matrix(as.vector(locs.unit), ncol=2)

matrix(Radial.basis(as.matrix(locs), center, 0.5), nrow=nrow(locs), byrow=TRUE)





x<- cbind( runif(100), runif(100))
center<- expand.grid( seq( 0,1,,5), seq(0,1,,5))
# coerce to matrix
center<- as.matrix(center)
PHI<- Radial.basis(x, center, delta=.5)







