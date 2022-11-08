library(argoFloats)
library(oce)
#> Loading required package: gsw
## 1. Get worldwide float-profile index, saving to ~/data/argo by default.
indexAll <- getIndex()
## 2. Narrow to a 30km-radius circle centred on Abaco Island, The Bahamas.
index <- subset(indexAll,
                circle=list(longitude=-77.06,latitude=26.54,radius=30))
#> Kept 41 cycles (0.00161%)
## 3. Get NetCDF files for these profiles, saving to ~/data/argo by default.
profiles  <- getProfiles(index)
## 4. Read the NetCDF files.
argos <- readProfiles(profiles)
#> Warning in readProfiles(profiles): Of 41 profiles read, 2 have >10% of conductivity values with QC flag of 4, signalling bad data.
#>     The indices of the bad profiles are as follows.
#>     3 8
#> Warning in readProfiles(profiles): Of 41 profiles read, 1 has >10% of pressure values with QC flag of 4, signalling bad data.
#>     The indices of the bad profiles are as follows.
#>     3
#> Warning in readProfiles(profiles): Of 41 profiles read, 4 have >10% of salinity values with QC flag of 4, signalling bad data.
#>     The indices of the bad profiles are as follows.
#>     3 6 7 13
#> Warning in readProfiles(profiles): Of 41 profiles read, 3 have >10% of temperature values with QC flag of 4, signalling bad data.
#>     The indices of the bad profiles are as follows.
#>     3 7 13
## 5. Examine QC flags, and set questionable data to NA.
argosClean <- applyQC(argos)
oldpar <- par(no.readonly=TRUE)
par(mfrow=c(1, 2))                     # want two-panel plot
par(mar=c(3.5, 2.0, 2.0, 2.0))         # maps do not get axis names
par(mgp=c(2,0.7,0))                    # tighten axes
## 6. Plot a map of profile locations.
plot(index, which="map", bathymetry=FALSE)
points(-77.06, 26.54, pch="*", cex=3)  # show centre of focus
mtext(paste(argosClean[["length"]], "profiles"), line=1.0)
## 7. Plot a TS diagram
par(mar=c(3.5, 3.5, 2.0, 1.0))         # increase left margin for name
plot(argosClean, which="TS")

#####

data("index")
plot(index, bathymetry=FALSE)
str(index)


# 6. function details
## 1. getIndex()
# getIndex(
#   filename = "core",
#   server = argoDefaultServer(),
#   destdir = argoDefaultDestdir(),
#   age = argoDefaultIndexAge(),
#   quiet = FALSE,
#   keep = FALSE,
#   debug = 0
# )
?getIndex
ai <- getIndex("core")
ai

## 2. subset()
# Subsetting by circle
aiCircle <- subset(ai, circle=list(longitude=-77.5, latitude=27.5, radius=50))
# Subsetting by polygon
lonPoly <- c(-76.5, -76.0, -75.5)
latPoly <- c(25.5, 26.5, 25.5)
aiPoly <- subset(ai, polygon=list(longitude=lonPoly, latitude=latPoly))
# Plotting the subsets together
CP <- merge(aiCircle, aiPoly)
plot(CP, bathymetry=FALSE)

??argoFloats::subset.
subset.
?argosubset
?points
plot(-80:-60, 20:40)
points(-76.5, 25.5, cex=2)
points(x=lonPoly, y=latPoly, cex=2)

######
