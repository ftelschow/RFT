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
# Developer notes:
# - Style guideline included
# - specify df for the different fields in the description
#------------------------------------------------------------------------------#

#' Computes the d-th EC density of different fields. Currently only "t" is implemented.
#'
#' @param x Numeric location x at which the EC density is to be evaluated
#' @param d Numeric d-th EC density
#' @param field String type of field for EC densities, currently only
#' "z" (Gaussian), "t" (Student-t) and "chi2" are supported
#' @param df Numeric degrees of freedom of the field.
#' @return value of EC density of the chosen field x.
#' @export
ECdensity <- function( x,
                       d,
                       field = "t",
                       df    = 1 ){
  if( field == "t" ){
    # Compute the EC densities for a t-field up to D = 3
    if( d == 0 ){
      1 - pt( x, df = df )

    }else if( d == 1 ){
      ( 2*pi )^(-1) * ( 1 + x^2 / df )^( -( df - 1 ) / 2 )

    }else if( d == 2 ){
      (2*pi)^(-3/2) * gamma( ( df + 1 ) / 2 ) / gamma( df / 2 ) / sqrt(df) *
        sqrt(2) * x * ( 1 + x^2 / df )^( -( df - 1 ) / 2 )

    }else if( d == 3 ){
      (2*pi)^(-2) * ( ( df - 1 ) / df * x^2 - 1 ) *
        ( 1 + x^2 / df )^( -( df - 1 ) / 2 )

    }else{
      stop( "Error: d must be smaller then 3. Higher dimensions
               are not yet implemented." )
    }

  }else if( field == "z" ){
    # Constant factor appearing in all EC densities
    constFunc = ( 2 * pi )^( -( d + 1 ) / 2 ) * exp( -x^2 / 2 )

    # Compute the EC densities for a z-field (gaussian) up to D = 3
    if( d == 0 ){
      1 - pnorm( x )

    }else if( d == 1 ){
      constFunc * 1

    }else if( d == 2 ){
      constFunc * 2 * x

    }else if( d == 3 ){
      constFunc * ( 4 * x^2 - 2 )

    }else{
      stop( "Error: d must be smaller then 4. Higher dimensions
               are not yet implemented." )
    }

  }else if( field == "chi2" ){
    # Compute the EC densities for a ch2-field up to D = 1
    if( d == 0 ){
      ifelse( x >= 0,
              pchisq( x, df = df, lower.tail = FALSE ),
              0 )

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
    stop( "Error: Input a valid field type." )
  }
}

#' Returns the expected Euler characterisitc curve for a specified field and its
#' Lipshitz-Killing-curvatures
#'
#' @param LKC Vector estimated LKCs of the UGRFs of the field.
#' @param field String type of field for EC densities, currently
#' "z" (Gaussian), "t" (Student-t) and "chi2" are supported. Default is "z".
#' @param df Numeric degrees of freedom of the field.
#' @return a function computing the EEC curve
#' @export
EEC <- function( LKC,
                 field  = "z",
                 df    = 0
                 ){

  #----- Approximate the tail distribution using the GKF and EECH
  if( length(LKC) == 2 ){
    # 1D case
    EECf <- function( u ){
              LKC[1] * EC_density( x    = u,
                                   d    = 0,
                                   field = field,
                                   df   = df ) +
              LKC[2] * EC_density( x    = u,
                                   d    = 1,
                                   field = field,
                                   df   = df )
    }

  }else if( length(LKC) == 3 ){
    # 2D case
    EECf <- function(u){
              LKC[1] * EC_density( x = u, d = 0, field = field, df = df ) +
              LKC[2] * EC_density( x = u, d = 1, field = field, df = df ) +
              LKC[3] * EC_density( x = u, d = 2, field = field, df = df )
    }

  }else if( length(LKC) == 4 ){
    # 3D case
    EECf <- function(u){
              LKC[1] * EC_density( x = u, d = 0, field = field, df = df ) +
              LKC[2] * EC_density( x = u, d = 1, field = field, df = df ) +
              LKC[3] * EC_density( x = u, d = 2, field = field, df = df ) +
              LKC[4] * EC_density( x = u, d = 3, field = field, df = df )
    }

  }else{
    stop( "Error: Not yet implemented for data with a domain of dimension
           greater than 3!" )
  }

  Vectorize( EECf )
}

#' Computes a threshold of the EEC curve. This can be used to control the FWER of
#' random fields.
#'
#' @param LKC Vector estimated LKCs of the UGRFs of the field.
#' @param alpha Numeric upper tail probabiltiy, which needs to be approximated.
#' @param field String type of field for EC densities, currently
#' "z" (Gaussian), "t" (Student-t) and "chi2" are supported. Default is "z".
#' @param df Numeric degrees of freedom of the field.
#' @param interval Vector with two components containing the upper and lower
#' bound for finding the root of the EEC curve. Default c(0,50).
#' @return Numeric the approximated 1-alpha quantile.
#' @export
EEC_threshold <- function( LKC,
                           alpha    = 0.05,
                           field    = "t",
                           df       = 0,
                           interval = c( 0, 50 )
                           ){
  # Get the EEC curve and subtract alpha
  EEC_function <- EEC( LKC, field = field, df )
  tailProb     <- function( u ){ EEC_function( u ) - alpha }

  # return the results
  list( threshold = uniroot( tailProb, interval )$root,
        EEC   = EEC_function,
        alpha = alpha,
        field = field
  )
}
