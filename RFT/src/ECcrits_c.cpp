/***********************************************
 *  This file contains c implementations of the
 *  Euler characteristic for different dimensions
 *  using the lower star/critical value trick
 *
***********************************************/
// Created by F. Telschow
// last modified: 5/2/2021

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

/* This is the C++ subroutine computing the change in Euler characterisitic in 1D
 * Input:
 *   - X: matrix containing as columns the observed field
 */
// [[Rcpp::export]]
int tupel2index( int x0, int x1, IntegerVector Dim ){
  // Get the domain of the input vector
  int D = Dim.length();

  // Create the output variable
  int out = x1 * Dim(0) + x0;
  
  return( out );
}

// [[Rcpp::export]]
int triple2index( int x0, int x1, int x2, IntegerVector Dim ){
  // Get the domain of the input vector
  int D = Dim.length();
  
  // Create the output variable
  int out = x2 * Dim(0) * Dim(1) + x1 * Dim(0) + x0;

  return( out );
}

/* This is the C++ subroutine computing the change in Euler characterisitic in 1D
 * Input:
 *   - X: matrix containing as columns the observed field
 */
// [[Rcpp::export]]
NumericMatrix ECcrit1D_C( NumericMatrix X ){
    // Get the number of columns of input, i.e., the number of realisations
    int K = X.nrow();
    int N = X.ncol();
    
    // Create a List with the number of columns of X
    NumericMatrix out( K, N );
    
    // Loop over the number if realisations
    for( int n = 0; n < N; n++ ){

        // Deal with the boundaries seperately
        if( X( 0, n ) > X( 1, n ) ){
            out( 0, n ) = -1;
        }
        
        if( X( K - 1, n ) > X( K - 2, n ) ){
            out( K-1, n ) = -1;
        }
        
        // Loop over the inner points of the domain
        for( int k = 1; k < K - 1; k++ ){
            // minima increase EC by one
            if( X( k-1, n ) > X( k, n ) && X( k, n ) < X( k+1, n ) ){
                out( k, n ) = 1;
            }
            
            // maxima decrease EC by one
            if( X( k-1, n ) < X( k, n ) && X( k, n ) > X( k+1, n ) ){
                out( k, n ) = -1;
            }
        }
    }

    return( out );
}


/* This is the C++ subroutine computing the change in Euler characterisitic in 1D
 * Input:
 *   - X: matrix containing as columns the observed field. The 2D Field is
 *        linearly stored.
 *   - dimX: the dimension of the original array
 */
// [[Rcpp::export]]
NumericVector ECcrit2D_C( NumericVector X, IntegerVector dimX ){
  // Get the number of columns of input, i.e., the number of realisations
  int K = dimX[0];
  int L = dimX[1];
  int N = dimX[2];
  
  // initialize the change in Euler characteristic as being 1, since 1 voxel is added
  int dec = 1;
  
  // Create a List with the number of columns of X
  NumericVector out( (K-2)*(L-2)*N );
  
  // Vector for indices
  IntegerVector ind(9);
  
  // Vector for neigbourhood
  double x[9];
  
  // Loop counter
  int i = 0;
  
  // Loop over the number if realisations
  for( int n = 0; n < N; n++ ){
    // Loop over the inner points of the domain
    for( int l = 1; l < L-1; l++ ){
      for( int k = 1; k < K-1; k++ ){
        /*
         * getting the indices of the 3x3 neighbourhood of
         * the voxel (ii,jj) and saving them in the ind_t array
         */
        ind[0] = triple2index( k-1, l-1, n, dimX );
        ind[1] = triple2index( k  , l-1, n, dimX );
        ind[2] = triple2index( k+1, l-1, n, dimX );
        ind[3] = triple2index( k-1, l,   n, dimX );
        ind[4] = triple2index( k  , l,   n, dimX );
        ind[5] = triple2index( k+1, l,   n, dimX );
        ind[6] = triple2index( k-1, l+1, n, dimX );
        ind[7] = triple2index( k  , l+1, n, dimX );
        ind[8] = triple2index( k+1, l+1, n, dimX );
        
        /* filling the x array using the indices of the
         * neighbourhood ind_t and the input array z
         */
        for( int ss = 0; ss < 9; ss++ )
        {
          x[ss] = X[ ind[ss] ];
        }
        
        /******************************************************************
         * compute the change in EC by checking which nD-faces are created
         ******************************************************************/
        // subtract number of edges
        IntegerVector ind1 = {1, 3, 5, 7};

        for( int ss = 0; ss < 4; ss++ ){
          // if new edge appears subtract 1 from dec
          if( x[ ind1[ss] ] > x[4] )
            dec -= 1;
        }
        
        // add number of faces
        if( x[0] > x[4] && x[1] > x[4] && x[3] > x[4] )
          dec += 1;
        if( x[1] > x[4] && x[2] > x[4] && x[5] > x[4] )
          dec += 1;
        if( x[3] > x[4] && x[6] > x[4] && x[7] > x[4] )
          dec += 1;
        if( x[5] > x[4] && x[7] > x[4] && x[8] > x[4] )
          dec += 1;
        
        // save the negative dec in the odd out locations (minus
        // sign convention used to recreate the ec curve)
        // Get the number of edges introduced by the center voxel
        out[i] = -dec;
        
        // reset dec to be 1
        dec = 1;
        // update index counter
        i += 1;
      }
    }
  }
  
  return( out );
}


/* This is the C++ subroutine computes the EC curves for given values
 * Input:
 *   - z: input array
 *   - cc: integer giving the connectivity (currently, only 8 is supported)
 *   - ADD all inputs with short discribtion
 */
// [[Rcpp::export]]
NumericMatrix ECcurve1D_C( NumericMatrix Y,
                           NumericMatrix dEC,
                           NumericVector u ){
    // Get the number of columns of input, i.e., the number of realisations
    int K = Y.nrow();
    int N = Y.ncol();
    int L = u.length();
    
    // Create the output matrix
    NumericMatrix out( L, N );
    
    // Loop over the number if realisations
    for( int n = 0; n < N; n++ ){
        // Get the indices where dEC is not equal zero
        LogicalVector ind_n = dEC( _ , n ) != 0;
        
        // Create the vector 1:K
        NumericVector lin(K); 
        for( int i = 0; i < K; i++ ){
            lin[i] = i;
        }
        
        // Get the rows where dEC is not equal zero
        NumericVector x_n = lin[ind_n];
            
        // Length of the vector x_n
        int M = x_n.length();
        
        // Loop over values in u
        for( int k = 0; k < L; k++ ){
            // Change ECu to zero
            NumericVector ECu = {0};
            
            // Check which critcal points are below u
            for( int m = 0; m < M; m++ ){
                // Add change in EC, if crit below u
                if( Y( x_n[m], n ) < u[k] ){
                    ECu[0] = ECu[0] + dEC( x_n[m], n );
                }
            }
            
            out( k, n ) = ECu[0];
        }
    }
    
    return( out );
}