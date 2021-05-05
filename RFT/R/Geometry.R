#------------------------------------------------------------------------------#
#
#     Functions for Geometric computations
#
#------------------------------------------------------------------------------#
# Contained functions:
#  - sphere_vol
#
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
sphere_vol <- function( j, d ){
  if( (d - 1 - j) %% 2 == 0 ){
    return( 2^{ j + 1 } * pi^{ j / 2 } * gamma( ( d + 1 ) / 2 ) /
                                 gamma( j + 1 ) / gamma( ( d + 1 - j ) / 2 ) )
  }else{
    return( 0 )
  }
}
