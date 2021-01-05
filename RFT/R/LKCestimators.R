################################################################################
####                                                                        ####
####             Selection of Estimators for the LKC of Gaussian fields     ####
####                                                                        ####
################################################################################
#
#  Implemented Estimators:
#     - direct LKC estimator by estimating the
#       Riemannian metric and computing the LKCs directly (1D and 2D)
#       as in Telschow Schwartzmann (2019)
#
################################################################# ###############
#' Computes the LKCs using the approximation of the defining integral method.
#' The integration is carried out using the package \code{polyCub}.
#'
#' @param R An array of dimension \code{c(length(x),length(y),n)} containing the
#'          realizations of the field.
#' @param coords a list containing as entries the x/y/z-Coordinates of the grid
#' on which the data is observed.
#' @param mask an object of type boundary as produced from the functions
#' \code{cont2bdry()} and \code{mask2bdry}.
#' @param h numeric. Offset for estimation of the derivatives.
#' @param plot logical. If true the integral domain is plotted.
#' @importFrom stats na.omit pnorm
#' @return a vector of length dim(R).
LKC_integral = function( R, coords = NULL, mask = NULL, h = 0.01, plot = F ){
  # Dimension of input residual array. Last coordinate is assumed to enumerate
  # obervations.
  dimR <- dim(R)
  D    <- length(dimR) - 1

  # Get default coordinate values
  if( is.null(coords) ){
    coords   <- list()
    coords$x <- 1:dimR[1]
    if( D > 1 ) coords$y <- 1:dimR[2]
    if( D > 2 ) coords$z <- 1:dimR[3]
  }else if( length( coords ) != D ){
    stop( "coords must include length(dim(R)) - 1 list entries." )
  }

  # number of samples
  nn = dimR[ D + 1 ]

  # get the coordinates as simple vectors
  x = coords$x
  if( D > 1 ){
    y = coords$y
  }
  if( D > 2 ){
    z = coords$z
  }

  # Switch estimation of LKCs by dimension
  if( D == 1 ){
    # Set the default mask, if no mask is provided
    if( is.null( mask ) ){
      mask <- list()
      mask[[1]] <- rbind( 1:length(x), x )
    }

    # Number of connected components
    nComp = length( mask )
    if( nComp == 0 ){
        return( function( u ) return( 0 ) )
    }

    #---------------------------------------------------------------------------
    # Compute L1
    #---------------------------------------------------------------------------
    # Loop over connected components
    L <- 0
    for( i in 1:nComp ){
      # observed loation on the connected component
      xx <- mask[[i]][2,]
      # get the derivative of the residuals
      dR <- apply( R[ mask[[i]][1,], ], 2,
                   FUN = function( yy ){
                             # Interpolate the function
                             fn <- stats::splinefun( x = xx,
                                                     y = yy,
                                                     method = "natural" )
                             # Obtain the derivative of the function
                             pracma::fderiv( f = fn,
                                             x = xx,
                                             n = 1,
                                             method = "central")
                         } )
      # get standard deviation of derivative
      dR.sd <- apply( dR , 1, stats::sd )

      # integrate the standard deviation of the derivative over the connected
      #  component using trapozoid rule
      L <- L + sum( diff(xx) / 2 * ( dR.sd[ 1:( length(dR.sd) - 1 ) ] + dR.sd[ 2:length(dR.sd) ] ) )
    }

  } else if( D == 2 ){
    # Initialize LKC vector
    L = rep( NA, D )
    names( L ) = c( "L1", "L2" )

    # Number of contour elements
    nCont = length( mask$bdry$bdry )
    if( nCont == 0 ){
        return( function( u ) return( 0 ) )
    }

    #---------------------------------------------------------------------------
    # Compute L1
    SC = vector( "list", nCont )

    # interpolate the residuals to the contour line from bdry object
    for( i in 1:nCont ){
      for( j in 1:nn ){
        if( j == 1 ){
          SC[[i]] = fields::interp.surface(
                      list( x = x, y = y, z = R[,,j] ),
                      cbind( mask$bdry$bdry[[i]]$x, mask$bdry$bdry[[i]]$y )
                    )
        }else{
          SC[[i]] = cbind( SC[[i]], fields::interp.surface(
            list( x = x, y = y, z = R[,,j] ),
            cbind( mask$bdry$bdry[[i]]$x, mask$bdry$bdry[[i]]$y ) ) )
        }
      }
    }

    # Gives the length of one component.
    L1 = function( X ){
      if( !is.matrix( X ) ) return( 0 ) # Pathological cases.
      if( nrow(X) == 1 ) return( 0 )
      sum( sqrt( rowSums( diff(X)^2 ) ) ) / sqrt( nn-1 )
    }

    # sums the length of all components wrt the Riemannian metric
    L[1] = sum( sapply( SC, L1 ) ) / 2

    #---------------------------------------------------------------------------
    # Compute L2
    volform <- function( ss ){
      N = dim( R )[3]
      dRx <- dRy <- matrix( NaN, dim(ss)[1], N )
      for( j in 1:N ){
        dRx[,j] <- ( fields::interp.surface(
          list( x = x, y = y, z = R[,,j] ),
          cbind( ss[,1] + h, ss[,2] )
        ) -
          fields::interp.surface(
            list( x = x, y = y, z = R[,,j] ),
            cbind( ss[,1] - h, ss[,2] )
          ) ) / 2 / h

        dRy[,j] <- ( fields::interp.surface(
          list( x = x, y = y, z = R[,,j] ),
          cbind( ss[,1], ss[,2] + h )
        ) -
          fields::interp.surface(
            list( x = x, y = y, z = R[,,j] ),
            cbind( ss[,1], ss[,2] - h )
          ) ) / 2 / h
      }

      sqrt( rowSums( dRx^2 ) * rowSums( dRy^2 ) - rowSums( dRx*dRy )^2 ) / ( N-1 )

    }

    # Use midpoint integration to get LKC2
    L[2] <- polyCub.midpoint( mask$bdry, f = volform, plot = TRUE )
  }

  return( L )
}
