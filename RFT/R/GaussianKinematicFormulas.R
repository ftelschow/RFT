#------------------------------------------------------------------------------#
#                                                                              #
#     Functions for Gaussian Kinematic formulas and thresholding               #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - EC_density
#  - EEC
#  - EEC_threshold
#
#------------------------------------------------------------------------------#
#' Computes the d-th EC density of different fields. Currently only "t" is implemented.
#'
#' @param x Numeric location x at which the EC density is to be evaluated
#' @param d Numeric d-th EX density
#' @param type String type of field for EC densities, currently only "t" (Student-t) and "z" (Gaussian) are supported
#' @param df Numeric degree of freedom of the field. For "t" this is a numeric.
#' @return value of 2d-EC density of a t-field of degree df at location x.
#' @export
EC_density <- function( x, d, type = "t", df = 1 ){
  # switch by type of the marginal distribution of the field
  if( type == "t" ){
    # compute the density depending on the dimension
    if( d == 0 ){
      1 - pt( x, df = df )
    }else if( d == 1 ){
      ( 2*pi )^(-1) * ( 1 + x^2 / df )^( -( df - 1 ) / 2 )
    }else if( d == 2 ){
      ( 2*pi )^( -3 / 2 ) * gamma( ( df + 1 ) / 2 ) / gamma( df / 2 ) /
      sqrt( df ) * sqrt(2) * x * ( 1 + x^2 / df )^( -( df - 1 ) / 2 )
    }else if( d == 3 ){
      ( 2*pi )^( -2 ) * ( ( df - 1 ) / df * x^2 - 1 ) *
      ( 1 + x^2 / df )^( -( df - 1 ) / 2 )
    }else{
      stop("Error: d must be smaller then 3. Higher dimensions are not yet implemented.")
    }

  }else if( type == "gauss" ){
    # function appearing in all densities
    constFunc = ( 2*pi )^( -( d + 1 ) / 2 ) * exp( -x^2 / 2 )

    # compute the density depending on the dimension
    if( d == 0 ){
      1 - pnorm( x )
    }else if( d == 1 ){
      constFunc * 1
    }else if( d == 2 ){
      constFunc * ( 2 * x )
    }else if( d == 3 ){
      constFunc * ( 4 * x^2 - 2 )
    }else{
      stop("Error: d must be smaller then 4. Higher dimensions are not yet implemented.")
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
#' @param LKC Vector estimated LKCs of the UGRFs of the field.
#' @param type String type of field.
#' @param df degrees of freedom.
#' @param alpha a constant to horizontally shift the EEC curve.
#' @return a function computing the EEC curve
#' @export
EEC <- function( LKC,
                 type  = "t",
                 df    = 0,
                 alpha = 0 ){
  ###### Approximate the tail distribution using the GKF and EECH
  if( length( LKC ) == 2 ){
    ### 1D case
    EECf <- function(u){
              LKC[1] * EC_density( x = u, d = 0,
                                   type = type, df = df ) +
              LKC[2] * EC_density( x = u, d = 1, type = type, df = df ) - alpha
    }

  }else if( length(LKC) == 3 ){
    ### 2D case
    EECf <- function(u){
              LKC[1] * EC_density( x = u, d = 0, type = type, df = df ) +
              LKC[2] * EC_density( x = u, d = 1, type = type, df = df ) +
              LKC[3] * EC_density( x = u, d = 2, type = type, df = df ) - alpha
    }

  }else if( length(LKC) == 4 ){
    ### 3D case
    EECf <- function(u){
              LKC[1] * EC_density( x = u, d = 0, type = type, df = df ) +
              LKC[2] * EC_density( x = u, d = 1, type = type, df = df ) +
              LKC[3] * EC_density( x = u, d = 2, type = type, df = df ) +
              LKC[4] * EC_density( x = u, d = 3, type = type, df = df ) - alpha
    }

  }else{
    stop( "Error: Not yet implemented for data with a domain of dimension greater than 3!")
  }

  Vectorize( EECf )
}

#' Computes a threshold of the EEC curve. This can be used to control the FWER of
#' random fields.
#'
#' @param alpha Numeric upper tail probabiltiy, which needs to be approximated.
#' @param LKC Vector estimated LKCs of the UGRFs of the field.
#' @param type String type of field.
#' @param df degrees of freedom.
#' @return Numeric the approximated 1-alpha quantile.
#' @export
EEC_threshold <- function( LKC,
                           alpha    = 0.05,
                           type     = "t",
                           df       = 0,
                           interval = c( 0, 50 )
                           ){

  tailProb <- EECf( LKC, type = df, df, alpha = alpha )

  list( threshold = uniroot( tailProb, interval )$root,
        EEC = Vectorize(function( u ) tailProb( u ) + alpha),
        alpha = alpha,
        type = type
  )
}
