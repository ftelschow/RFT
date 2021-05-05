#------------------------------------------------------------------------------#
#                                                                              #
#     Script demonstrating the features of the package                         #
#                                                                              #
#------------------------------------------------------------------------------#
# Chapter 1: LKC_integral
# Chapter 2: Getting the threshold from the GKF
#  
#
#------------------------------------------------------------------------------#
.rs.restartR()
# clean workspace
rm( list = ls() )

# required packages
require( SCBfda )
require( spatstat )
require( polyCub )
require( LKC )

#------------------------------------------------------------------------------#
# Chapter 1: LKC estimation
#------------------------------------------------------------------------------#
#  Generate the data
x = seq( 1, 70, by = 1 )
coords <- list( x = x, y = x )
test <- SquaredExp2DNoise( 100, x = x )

# Generate the mask
mask = matrix( 0, length(x), length(x) )
mask[ 11:60, 11:60 ] <- 1

# Get Euler Characteristic of the mask/domain. This is LKC 0.
L0 <- EulerChar( mask, 0.5, connectivity = 8 )
  
# Get the boundary from the mask
cont <- mask2bdry( coords, mask )

# Plot mask and Contour lines
image( x, x, mask)
contour( list( x = x,
               y = x, z = mask ),
         levels = 1, nlevels = 1, add = TRUE )

# Get the normalized residuals
tmp = apply( test, 1:2, scale )
R   = aperm( tmp, c( 2, 3, 1 ) )

# Estimate the LKCs
Ltest = LKC_integral( coords, cont, R )
Ltest

#-------------------------------------------------------------------------------
# Example with complicated mask
#-------------------------------------------------------------------------------
# Generate the mask
mask = matrix( 0, length(x), length(x) )
mask[ 11:60, 11:60 ] <- 1
mask[ 30:40, 30:40 ] <- 0
mask[ 15:20, 15:20 ] <- 0
mask[ 35:37, 35:37 ] <- 1

# Get the boundary from the mask
bdry <- mask2bdry( coords, mask )

plot( bdry$mask, col = 2)
plot( bdry$bdry, add = T)

# Get second method for the boundary from the mask
cont = contourLines( list( x = x,
                           y = x, z = mask ),
                     levels = 1, nlevels = 1 )
bdry2 <- cont2bdry( coords, cont, hole = bdry$hole )

image( x = x, y = x, mask )
plot( bdry2$bdry, add = T)

# Get Euler Characteristic of the mask/domain. This is LKC 0.
L0 <- EulerChar( mask, 0.5, connectivity = 8 )

# Estimate the LKCs
LKC_integral( coords, bdry, R )
LKC_integral( coords, bdry2, R ) # There seems to be still a bug with the holes.
                               # it works for simply connected domains though
Ltest

tmp = bdry2
tmp <- reverse_contour( tmp, 3 )
LKC_integral( coords, tmp, R, plot = T )

##### Bug fixing
# this seems to be alright. Not sure why the integral does not see the latter
# two parts as holes.
sapply( bdry$bdry$bdry, find_holes ) 
sapply( bdry2$bdry$bdry, find_holes )
sapply( tmp$bdry$bdry, find_holes )

#------------------------------------------------------------------------------#
# Chapter 2: GKF threshold
#------------------------------------------------------------------------------#
