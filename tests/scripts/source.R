install_github("ArgoCanada/argoFloats", ref="develop")

library(oce)
library(ocedata)
library(argoFloats)
data(index)
str(index)
plot(index, bathymetry=TRUE)
install.packages("ncdf4")




ai <- getIndex("core")

lonPoly <- c(-76.5, -76.0, -75.5)
latPoly <- c(25.5, 26.5, 25.5)
aiPoly <- subset(ai, polygon=list(longitude=lonPoly, latitude=latPoly))

plot(aiPoly, bathymetry=FALSE)
from <- as.POSIXct("2020-09-23", tz="UTC")
to <- as.POSIXct("2020-10-25", tz="UTC")

index2 = subset(aiPoly, time = list(from=from, to=to))

plot(index2, bathymetry=FALSE,  asp=1/cos(mean(range(unlist(index2[["latitude"]]), na.rm=TRUE))*pi/180),  mgp=getOption("oceMgp")
)

#
lonlim <- c(-70, -64,-10)
latlim <-c(40,35,35)
index1 <- subset(ai, section=list(longitude=lonlim, latitude=latlim, width=100))
#subset by time
from <- as.POSIXct("2020-09-23", tz="UTC")
to <- as.POSIXct("2020-10-25", tz="UTC")
index2 <- subset(index1, time=list(from=from, to=to))
plot(index2, bathymetry=FALSE,  asp=1/cos(mean(range(unlist(index2[["latitude"]]), na.rm=TRUE))*pi/180),  mgp=getOption("oceMgp")
)
points(lonlim, latlim, pch=21, col="black", bg="red", type="o")
plot(section, which="map", col="tan")
par(oldpar)

prof2 <- getProfiles(index2)

argos <- readProfiles(prof2)
?applyQC

a <- c(1i, 1); b <- c(-1i, 1)

T <- cbind(c(2,3), c(-3, 2))
T %*% c(1,2)

eigen(T)




