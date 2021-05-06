#------------------------------------------------------------------------------#
#                                                                              #
#     Functions for the Gaussian Kinematic formulas and thresholding           #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - ECdensity
#  - GKFthreshold
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline included
# - specify df for the different fields in the description
#------------------------------------------------------------------------------#

#' Computes the d-th EC density of different fields. Currently only "t" is implemented.
#'
#' @param x Numeric location x at which the EC density is to be evaluated
#' @param d Numeric d-th EX density
#' @param type String type of field for EC densities, currently only "t" (Student-t) and "z" (Gaussian) are supported
#' @param df Numeric degree of freedom of t-field
#' @return value of 2d-EC density of a t-field of degree df at location x.
#' @export
ECdensity <- function( x, d, type = "t", df = 1 ){
  if( type == "t" ){
    if( d == 0 ){
      1 - pt( x, df = df )
    }else if( d == 1 ){
      (2*pi)^(-1) * ( 1 + x^2/df )^( -(df-1)/2 )
    }else if( d == 2 ){
      (2 * pi)^(-3/2) * gamma( (df + 1) / 2 ) / gamma( df/2 ) / sqrt(df) * sqrt(2) * x * ( 1 + x^2/df )^( -(df-1)/2 )
    }else if( d == 3 ){
      (2*pi)^(-2) * ( (df-1)/df*x^2 - 1 ) * ( 1 + x^2/df )^( -(df-1)/2 )
    }else{
      stop("Error: d must be smaller then 4. Higher dimensions are not yet implemented.")
    }
  }else if( type == "gauss" ){
    constFunc = (2*pi)^( -( d + 1 ) / 2 ) * exp( -x^2 / 2 )
    if( d == 0 ){
      1 - pnorm( x )
    }else if( d == 1 ){
      constFunc * 1
    }else if( d == 2 ){
      constFunc * (2*x)
    }else if( d == 3 ){
      constFunc * ( 4 * x^2 - 2 )
    }else{
      stop("Error: d must be smaller then 4. Higher dimensions are not yet implemented.")
    }
  }else if( field == "T" ){
    # Compute the EC densities for Hotelling T-field for vectors up to D = 3.
    if( df[1] <= 3 ){
      out = 0
      if( x > 0 ){
        for( j in 0:( df[1] - 1 ) ){
          out = out + sphere_vol( j, df[1] ) * EC_density( sqrt( x ),
                                                           d + j,
                                                           field = "t",
                                                           df    = df[2] )
        }
      }else{
        out = 1
      }
      
      return( out )
    }else{
      stop("Error: d must be smaller then 4. Higher dimensions are not yet implemented.")
    }
  }else if( field == "chi2" ){
    # Compute the EC densities for a ch2-field up to D = 1
    if( d == 0 ){
      ifelse( x >= 0,
              pchisq( x, df = df, lower.tail = FALSE ),
              1 )

    }else if( d == 1 ){
      ifelse( x >= 0,
              x^( ( df - 1 ) / 2 ) * exp( -x / 2 ) / sqrt( 2 * pi ) /
                2^( df / 2 - 1 ) / gamma( df / 2 ),
              0 )

    }else{
      stop( "Error: d must be smaller then 4. Higher dimensions are
             not yet implemented." )
    }

  }else{
    stop("Error: Input a valid field type.")
  }
}

#' Approximates the upper tail probabilities of different fields using the corresponding Gaussian
#' kinematic formulas.
#' It approximates the tail probability of the absolute value by twice the
#' usual tail probability of the process.
#'
#' @param alpha Numeric upper tail probabiltiy, which needs to be approximated.
#' @param LKC Vector estimated LKCs of the UGRFs of the field.
#' @param type String type of field.
#' @param df degrees of freedom.
#' @return Numeric the approximated 1-alpha quantile.
#' @export
GKFthreshold <- function( alpha = 0.05, LKC, type = "z", df = 0, interval = c(0, 100) ){
  ###### Approximate the tail distribution using the GKF and EECH
  if( length( LKC ) == 2 ){
    ### 1D case
    tailProb <- function(u){
      LKC[1] * ECdensity( x = u, d = 0,
                          type = type, df = df ) +
      LKC[2] * ECdensity( x = u, d = 1, type = type, df = df ) - alpha
    }
  }else if( length(LKC) == 3 ){
    ### 2D case
    tailProb <- function(u){
        LKC[1] * ECdensity( x = u, d = 0, type = type, df = df ) +
        LKC[2] * ECdensity( x = u, d = 1, type = type, df = df ) +
        LKC[3] * ECdensity( x = u, d = 2, type = type, df = df ) - alpha
    }
  }else if( length(LKC) == 4 ){
    ### 3D case
    tailProb <- function(u){
      LKC[1] * ECdensity( x = u, d = 0, type = type, df = df ) +
      LKC[2] * ECdensity( x = u, d = 1, type = type, df = df ) +
      LKC[3] * ECdensity( x = u, d = 2, type = type, df = df ) +
      LKC[4] * ECdensity( x = u, d = 3, type = type, df = df )
      - alpha
    }
  }else{
    stop( "Error: Not yet implemented for data with a domain of dimension greater than 3!")
  }
  list( threshold = uniroot( tailProb, interval = interval )$root,
        EEC       = Vectorize(function(u) tailProb(u) + alpha),
        alpha     = alpha,
        type      = type )
}
