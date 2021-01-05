#------------------------------------------------------------------------------#
#                                                                              #
#     Miscalaneous functions for masking and creating boundaries               #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - mask2bndry
#  - cont2bdry
#  - find_holes
#  - reverse_contour
#
#------------------------------------------------------------------------------#
#' This function takes as input a binary mask and returns a boundary object as
#' required for the LKC estimation methods.
#'
#' @param x x-Coordinates of the grid on which the mask is defined.
#' @param y y-Coordinates of the grid on which the mask is defined.
#' @param mask the mask object. Must be binary, i.e. only contain 0's and 1's.
#' @return an object of class contour.
mask2bdry = function( mask, coords, n = 20 ){
    x = coords$x
    y = coords$y
    # create an owin object from the mask
    mask  <- owin( xrange = range(x),
                   yrange = range(y),
                   mask   = t(matrix( as.logical(mask), dim(mask)[1], dim(mask)[2] )) )
    bdry <- as.polygonal( mask )

    hole <- rep( FALSE, length( bdry$bdry ) )

    for( j in 1:length( bdry$bdry ) ){
      # Figure out whether the line is a hole or not
      tmp = find_holes( bdry$bdry[[j]] )
      hole[j] = tmp$hole
      xvals   = tmp$xvals
      yvals   = tmp$yvals

      # Linear interpolation to increase resolution
      lambda <- seq( 0, 1, length.out = n )
      xvals_n <- yvals_n <- NULL
      for( k in 1:(length( xvals ) - 1 ) ) {
        # Get the linear interpolation
        tmp = rbind( (1-lambda) * xvals[k], (1-lambda) * yvals[k] ) +
              rbind( lambda * xvals[k+1], lambda * yvals[k+1] )
        xvals_n <- c( xvals_n, tmp[1,] )
        yvals_n <- c( yvals_n, tmp[2,] )
      }
      bdry$bdry[[j]]$x <- xvals_n
      bdry$bdry[[j]]$y <- yvals_n
  }

  list( mask = mask, bdry = bdry, hole = hole )
}

#------------------------------------------------------------------------------#
#' This function takes as input the output of the \code{contourLines()} function
#' and creates an owin based data object for LKC estimation from it.
#' Note that the indication for holes needs to be done manually. By default all
#' contour lines are interpreted as filled areas!
#'
#' @param x x-Coordinates of the grid on which the mask is defined.
#' @param y y-Coordinates of the grid on which the mask is defined.
#' @param cont output from \code{contourLines()}.
#' @return an object of class contour.
cont2bdry = function( coords, cont, hole ){
  x = coords$x
  y = coords$y
  # number of contourlines in cont
  nCont <- length(cont)
  # initialize output object
  out = list()
  out$bdry <- 1
  out$hole <- hole

  # create an poly feature owin object from the cont
  for( j in 1:nCont ){
    if( j == 1 ){
      if( out$hole[j] != find_holes( cont[[j]] )$hole ){
        poly <- cbind( rev( cont[[j]]$x[-1] ), rev( cont[[j]]$y[-1] ) )
      }else{
        poly <- cbind( cont[[j]]$x[-1], cont[[j]]$y[-1] )
      }
      out$bdry <- owin( x = range(x), y = range(y), poly = poly )
    }else{
      out$bdry$bdry[[j]] <- list()
      if( out$hole[j] != find_holes( cont[[j]] )$hole ){
        out$bdry$bdry[[j]]$x <- rev( cont[[j]]$x[-1] )
        out$bdry$bdry[[j]]$y <- rev( cont[[j]]$y[-1] )
      }else{
        out$bdry$bdry[[j]]$x <- cont[[j]]$x[-1]
        out$bdry$bdry[[j]]$y <- cont[[j]]$y[-1]
      }

    }
  }

  return( out )
}

#------------------------------------------------------------------------------#
#' This function takes as input the output of the \code{contourLines()} function
#' and creates an owin based data object for LKC estimation from it.
#' Note that the indication for holes needs to be done manually. By default all
#' contour lines are interpreted as filled areas!
#'
#' @param bdry output from \code{contourLines()}.
#' @return an object of class contour.
find_holes = function( bdry ){
  # initialize output logical vector for holes
  hole  <- TRUE

    xvals <- bdry$x
    xvals <- c( xvals, xvals[1] )
    yvals <- bdry$y
    yvals <- c( yvals, yvals[1] )
    if( sum( diff(xvals) * ( yvals[-length(yvals)]+yvals[-1] ) ) < 0 ){
      hole <- FALSE
    }
  list( xvals = xvals, yvals = yvals, hole = hole )
}

#------------------------------------------------------------------------------#
#' This function reverses the order of specified contours in a bdry object. This
#' is required, since  the output of the \code{contourLines()} function does not
#' necessarily has the coorect orientation to indicate that the contour belongs
#' to a hole.
#'
#' @param bdry output from \code{contourLines()}.
#' @return an object of class contour.
reverse_contour = function( bdry, indices ){
  for( j in indices ){
    bdry$bdry$bdry[[j]]$x <- rev( bdry$bdry$bdry[[j]]$x )
    bdry$bdry$bdry[[j]]$y <- rev( bdry$bdry$bdry[[j]]$y )

    bdry$hole[j] <- !bdry$hole[j]
  }

  return( bdry )
}
