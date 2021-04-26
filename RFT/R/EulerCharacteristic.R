#------------------------------------------------------------------------------#
#
#     Functions for Euler characteristic computations
#
#------------------------------------------------------------------------------#
# Contained functions:
#  - EulerChar
#
#------------------------------------------------------------------------------#

#' Estimates the Euler characteristic of the exceedance set A(u)={ s in S | f(s) >u }
#' of a function f:S->R. 4 or 8 connectivity  is currently used and it is
#' implemented up to dimension 2.
#'
#' @param f Vector/Matrix observation of a random function on a grid
#' @param u Vector exceedance values to be computed
#' @param connectivity Integer amount of points considered as neighbors
#' @return Vector of Euler characteristics of the exceedance sets for the values u
#' @export
EulerChar <- function( f, u, connectivity = 8 ){
  # Make vector input a matrix
  if( is.vector( f ) ){
    f <- matrix( f, 1, length( f ) )
  }

  sz = dim( f )
  N  = ifelse( sz[1] == 1, 1, length( sz ) )

  EC <- rep( 0, length(u) )

  for( j in length(u):1 ){
    A = ( f >= u[j] )

    if( N == 1 ){
      vertices <- sum( A )
      edges    <- sum( A[ 1:sz[2] - 1 ] & A[ 2:sz[2] ] )
      EC[j]    <- vertices - edges

    }else if( N == 2 ){
      if( connectivity == 4 ){
        vertices <- sum( A )
        edges    <- sum( A[1:(sz[1]-1),] & A[2:sz[1],] ) + sum( A[,1:(sz[2]-1)] & A[,2:sz[2]] )
        faces    <- sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:(sz[1]),1:(sz[2]-1)] & A[1:(sz[1]-1),2:sz[2]] & A[2:sz[1],2:sz[2]] )
        EC[j]    <- vertices - edges + faces ;

      }else if( connectivity == 8 ){
        vertices <- sum( A )
        edges    <- sum( A[1:(sz[1]-1), ] & A[2:sz[1], ] ) + sum( A[ ,2:sz[2]] & A[ ,1:(sz[2]-1)] ) +
          sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:sz[1],2:sz[2]] ) + sum( A[2:sz[1],1:(sz[2]-1)] & A[1:(sz[1]-1),2:sz[2]] )
        # Summands left to right correspond to the following faces
        #   o--o    o  x   o--o    x  o
        #   \ /     \ \      \\     / \
        #   o  x    o--o   x  o    o--o
        p.faces  <- sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:sz[1],1:(sz[2]-1)] &  A[1:(sz[1]-1),2:sz[2]] ) +
          sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:sz[1],1:(sz[2]-1)] &  A[2:sz[1],2:sz[2]] ) +
          sum( A[1:(sz[1]-1),2:sz[2]]     & A[2:sz[1],2:sz[2]]     &  A[1:(sz[1]-1),1:(sz[2]-1)] ) +
          sum( A[1:(sz[1]-1),2:sz[2]]     & A[2:sz[1],2:sz[2]]     &  A[2:sz[1],1:(sz[2]-1)] )
        # Note that we only count 4 instead of 5 verices and 6 instead of 8 edges. Thus, we need to subtract a face if the following scheme appears
        # o---o
        # \ x \
        # o---o
        n.faces  <- sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:(sz[1]),1:(sz[2]-1)] & A[1:(sz[1]-1),2:sz[2]] & A[2:sz[1],2:sz[2]] )
        EC[j]    <- vertices - edges + p.faces - n.faces ;

      }else{
        stop("For 2D fields 'connectivity' has to have the value 4 or 8.")
      }

    }else{
      stop("N >= 2 currently not implemented")
    }
  }
  return( EC )
}


#' Estimates the Euler characteristic of the exceedance set A(u)={ s in S | f(s) >u }
#' of a function f:S->R. 4 or 8 connectivity  is currently used and it is
#' implemented up to dimension 2.
#'
#' @param f Vector/Matrix observation of a random function on a grid
#' @param D Integer dimension of the domain
#' @return A function computing the Euler characteristic curve
#' @export
EulerCharC <- function( f, Dim = NULL ){
  # Make vector input a matrix
  if( class( f ) != "RandomField"  ){
    break( "Error: the input must be of class RandomField." )
  }
  
  # Get the EC of the underlying domain
  L0 = 1
  
  # Get the dimension of the domain
  D = f$D
  
  if( D == 1 ){
    dEC = ECcrit1D_C( f$values )
    
    EC = function( u ) L0 + ECcurve1D_C( f$values, dEC, u )
    
  }else if( D == 2 ){
    tildeA = array( -Inf, dim = Dim + c( 2, 2, 0) )
    tildeA[ (1:Dim[1]) + 1,  (1:Dim[2]) + 1, ] = array( c(f$values), dim = Dim )
    dEC = ECcrit2D_C( tildeA, dim(tildeA) )

    EC  = function( u ) L0 + ECcurve1D_C( matrix( f$values, f$nloc, f$N ),
                                          matrix( dEC, f$nloc, f$N ),
                                          u )
    dEC = array( c(dEC), dim = Dim )
    
  }else if( D == 3 ){
    
  }else{
    break( "Error: dimension of the domain of the RandomField must be smaller
           than 4." )
  }

  # Create the output list
  out <- list()
  out[["EC"]]  <- EC
  out[["dEC"]] <- dEC
  
  return( out )
}