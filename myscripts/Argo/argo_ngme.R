dat <- read.table('myscripts/Argo/ocean_data.txt',header=TRUE)
str(dat)

max.edge    = 1
bound.outer = 5
locY = unique(cbind(dat$lon, dat$lat))
mesh = inla.mesh.2d(loc=locY,
                    # the inner edge and outer edge
                    max.edge = c(1,5)*max.edge,
                    # offset extension distance inner and outer extenstion
                    offset = c(max.edge, bound.outer)
)
plot(mesh)
points(locY, col = "red")

fem_mesh = inla.mesh.fem(mesh)
Ce <- fem_mesh$c1 #<phi_i, phi_j>
C <- fem_mesh$c0 #approximation of Ce
G <- fem_mesh$g1
A <- inla.spde.make.A(mesh, locY) #dim(A) = data loc * vertices

#There are some columns in the projector matrix all of whose elements equal zero:
table(colSums(A) > 0)
#These columns correspond to triangles with no point location inside.
#These columns can be dropped.
P <-mesh$loc
FV<-mesh$graph$tv

