



#include <R.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>



void crispify( double *A_cc, double *A_cs, double *A_s, double *A_z, 
              double *B_cc, double *B_cs, double *B_s, double *B_z, 
              double *Z_in, double *Z_out, double *alpha, double *flatten,
              int *start_self_ref, int *len_ssr, 
			  int *z_lags_vec, int *z_rep_in, int *n_Zd,
			  int *d_A, int *d_B, int *t_tau, int *stnd_d, int *tss, int *geodesic, double *logTen_cutoff ) {
	/// output MX will be d_a x d_b -- i.e., A is treated is X, B as Y in covariance Mx
	int i, ii, j, n = d_A[0], ktn = d_B[0] ;
    int tau = t_tau[0] ;
    int i_x_tau ;
    double wflatten = flatten[0] ;
    double malpha = -alpha[0] ;
	double term1, term2 ;
    double w_vec[ ktn ] ;
	double n_vec_prod ;
    double BIG_SUM_1, BIG_SUM_2 ;
    double BIG_SUM_VEC[ ktn ] ;
    //double inv_dstnd ; // commented out 20191204 -- v 1.01
    //double ktn_x_inv_dstnd ; // commented out 20191204 -- v 1.01
    double inv_wstnd ; 
    double inv_wstnd_x_wflatten ;
    double sparse_cutoff ;
    sparse_cutoff = pow( 10, logTen_cutoff[0] ) ;
    
    //printf( "%d %d \n" , n, ktn );
    
    int sparse_ndx_vec[ ktn ] ;
    int snv_k ;
	
	//int size_Zd = tau * n_Zd[0] ; // commented out 20191204 -- v 1.01
	int this_ndx ;
	int z_in_key ;
	
	int tau_zin_key_vec[ ktn ] ;
	for( j=0; j<ktn; j++ ) {
		tau_zin_key_vec[j] = tau * z_rep_in[ j ] + z_lags_vec[ j ] ;
	}
    
    for( j=0; j<ktn; j++ ) { BIG_SUM_VEC[j] = 1 ; } ////////////////////////////////////////////////
    
	if( stnd_d[0] == 1 ) {
		
		if( geodesic[0] == 1 ) {
			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				for( i=0; i<n; i++ ) { ////// for each site
					n_vec_prod = A_cc[i]*B_cc[j] + A_cs[i]*B_cs[j] + A_s[i]*B_s[j] ;
					if(  n_vec_prod < -1 ) { n_vec_prod = -1 ; }
					if(  n_vec_prod > 1 ) { n_vec_prod = 1 ; }
					//Rprintf( "%g \t" , n_vec_prod );
					term1 = 6378.1 * acos( n_vec_prod ) ;
					
					//term2 = pow(  pow( term1, 2 )  +   pow( A_z[i] - B_z[j], 2 ) , 0.5 ) ; // DZ edit point --- otherwise sites seperated by lags will b over amplified !!
					term2 = term1 ;
					
					BIG_SUM_VEC[j] = BIG_SUM_VEC[j] + term2 ;
					//Rprintf( "%g \t" , term1 );
				}
				//Rprintf( "\n" );
			}
		
		} else { ////// not geodesic

			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				for( i=0; i<n; i++ ) { ////// for each site
					term1 = pow(  pow( A_cc[i] - B_cc[j], 2 ) + pow( A_cs[i] - B_cs[j], 2 ),   0.5 ) ;
					BIG_SUM_VEC[j] = BIG_SUM_VEC[j] + term1 ;
					//Rprintf( "%g \t" , term1 );
				}
				//Rprintf( "\n" );
			}
			
		} // end geodesic or not
		
		for( j=0; j<ktn; j++ ) { BIG_SUM_VEC[j] = pow( BIG_SUM_VEC[j], -1 ) * n ; }
	
	} // end stnd_d ---- if NOT standardize D, BIG_SUM_VEC already set to 1, so no else necessary
    
	
	
    //for( j=0; j<ktn; j++ ) { Rprintf( "%g \t" , BIG_SUM_VEC[j] ) ; }
    
	for( i=0; i<n; i++ ) { ////// for each site
		
        i_x_tau = i * tau ;
		
		if( geodesic[0] == 1 ) {
			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				n_vec_prod = A_cc[i]*B_cc[j] + A_cs[i]*B_cs[j] + A_s[i]*B_s[j] ;
				if(  n_vec_prod < -1 ) { n_vec_prod = -1 ; }
				if(  n_vec_prod > 1 ) { n_vec_prod = 1 ; }
				term1 = 6378.1 * acos( n_vec_prod ) ;
				term2 = pow(  pow( term1, 2 )  +   pow( A_z[i] - B_z[j], 2 ) , 0.5 ) ;
				w_vec[ j ] = term2 * BIG_SUM_VEC[j] ;
				//Rprintf( "%g \t" , w_vec[ j ] );
				//Rprintf( "%g \t" , term1 );
			}
			
		} else {
			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				term1 = pow(  pow( A_cc[i] - B_cc[j], 2 )  +  pow( A_cs[i] - B_cs[j], 2 )  +  pow( A_z[i] - B_z[j], 2 ),   0.5 ) ;
				term2 = term1 ;
				w_vec[ j ] = term2 * BIG_SUM_VEC[j] ;
				//Rprintf( "%g \t" , w_vec[ j ] );
				//Rprintf( "%g \t" , term1 );
			}
			
				
		}
		
		
        //Rprintf( "\n" );
		
        for( j=0; j<ktn; j++ ) {
            w_vec[ j ] = exp( malpha * w_vec[ j ] ) ;
            //Rprintf( "%g \t" , w_vec[ j ] );
        }
        //Rprintf( "\n" );
        
        if( start_self_ref[0] > -1 ) {
            for( j=0; j<len_ssr[0]; j++ ) { ///// start_self_ref must be zero-base vector
                w_vec[ start_self_ref[j]+i ] = 0 ;
                //Rprintf( "%g \t" , w_vec[ j ] );
            }
        }
        //Rprintf( "\n" );
        
        
        BIG_SUM_1 = 0 ;
        for( j=0; j<ktn; j++ ) {
            BIG_SUM_1 = BIG_SUM_1 + w_vec[ j ] ;
        }
        //Rprintf( "%g \t" , BIG_SUM_1 );
        
        if( BIG_SUM_1 == 0 ) { inv_wstnd = 1 ; } else { inv_wstnd = pow( BIG_SUM_1, -1 ) ; }
        //inv_wstnd = pow( BIG_SUM_1, -1 ) ;
		//Rprintf( "%g \t" , BIG_SUM_1 );
        
        inv_wstnd_x_wflatten = inv_wstnd * wflatten ;
        for( j=0; j<ktn; j++ ) {
            w_vec[ j ] = w_vec[ j ] * inv_wstnd_x_wflatten ;
            //Rprintf( "%g \t" , w_vec[ j ] );
        }
        //Rprintf( "\n" );
        
        snv_k = 0 ;
        for( j=0; j<ktn; j++ ) {
            if( !isnan(w_vec[ j ]) ) {
                if( w_vec[ j ] > sparse_cutoff ) {
                    sparse_ndx_vec[ snv_k ] = j ;
                    snv_k = snv_k + 1 ;
					//Rprintf( "%g \t" , w_vec[ j ] );
                }
				//Rprintf( "\n" );
            }
			
		}
        
        
        for( ii=tss[0] ; ii<tss[1] ; ii++ ) {
            BIG_SUM_2 = 0 ;
            for( j=0; j<snv_k; j++ ) {
				this_ndx = sparse_ndx_vec[j] ; 
				
				//z_in_key = tau * z_rep_in[ this_ndx ] + ( z_lags_vec[ this_ndx ] + ii ) ;
				
				z_in_key = tau_zin_key_vec[ this_ndx ]  + ii ;
				
				//if( z_in_key >= 0 & z_in_key <= size_Zd ) {
					BIG_SUM_2 = BIG_SUM_2 + w_vec[ this_ndx ] * Z_in[ z_in_key ] ; /// can move flatten here
				//}
            }
            Z_out[ ii + i_x_tau ] = BIG_SUM_2 ;
        }
		
	}
}






void crispify_SLOW( double *A_cc, double *A_cs, double *A_s, double *A_z, 
              double *B_cc, double *B_cs, double *B_s, double *B_z, 
              double *Z_in, double *Z_out, double *alpha, double *flatten,
              int *start_self_ref, int *len_ssr, 
			  int *z_lags_vec, int *z_rep_in, int *n_Zd,
			  int *d_A, int *d_B, int *t_tau, int *stnd_d, int *tss, int *geodesic, double *logTen_cutoff ) {
	/// output MX will be d_a x d_b -- i.e., A is treated is X, B as Y in covariance Mx
	int i, ii, j, n = d_A[0], ktn = d_B[0] ;
    int tau = t_tau[0] ;
    int i_x_tau ;
    double wflatten = flatten[0] ;
    double malpha = -alpha[0] ;
	double term1, term2 ;
    double w_vec[ ktn ] ;
	double n_vec_prod ;
    double BIG_SUM_1, BIG_SUM_2 ;
    double BIG_SUM_VEC[ ktn ] ;
    //double inv_dstnd ; // commented out 20191204 -- v 1.01
    //double ktn_x_inv_dstnd ; // commented out 20191204 -- v 1.01
    double inv_wstnd ; 
    double inv_wstnd_x_wflatten ;
    double sparse_cutoff ;
    sparse_cutoff = pow( 10, logTen_cutoff[0] ) ;
    
    //printf( "%d %d \n" , n, ktn );
    
    int sparse_ndx_vec[ ktn ] ;
    int snv_k ;
	
	int size_Zd = tau * n_Zd[0] ;
	int this_ndx ;
	int z_in_key ;
    
    for( j=0; j<ktn; j++ ) { BIG_SUM_VEC[j] = 1 ; } ////////////////////////////////////////////////
    
	if( stnd_d[0] == 1 ) {
		
		if( geodesic[0] == 1 ) {
			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				for( i=0; i<n; i++ ) { ////// for each site
					n_vec_prod = A_cc[i]*B_cc[j] + A_cs[i]*B_cs[j] + A_s[i]*B_s[j] ;
					if(  n_vec_prod < -1 ) { n_vec_prod = -1 ; }
					if(  n_vec_prod > 1 ) { n_vec_prod = 1 ; }
					//Rprintf( "%g \t" , n_vec_prod );
					term1 = 6378.1 * acos( n_vec_prod ) ;
					
					//term2 = pow(  pow( term1, 2 )  +   pow( A_z[i] - B_z[j], 2 ) , 0.5 ) ; // DZ edit point --- otherwise sites seperated by lags will b over amplified !!
					term2 = term1 ;
					
					BIG_SUM_VEC[j] = BIG_SUM_VEC[j] + term2 ;
					//Rprintf( "%g \t" , term1 );
				}
				//Rprintf( "\n" );
			}
			
		} else { ////// not geodesic
			
			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				for( i=0; i<n; i++ ) { ////// for each site
					term1 = pow(  pow( A_cc[i] - B_cc[j], 2 ) + pow( A_cs[i] - B_cs[j], 2 ),   0.5 ) ;
					BIG_SUM_VEC[j] = BIG_SUM_VEC[j] + term1 ;
					//Rprintf( "%g \t" , term1 );
				}
				//Rprintf( "\n" );
			}
			
		} // end geodesic or not
		
		for( j=0; j<ktn; j++ ) { BIG_SUM_VEC[j] = pow( BIG_SUM_VEC[j], -1 ) * n ; }
		
	} // end stnd_d ---- if NOT standardize D, BIG_SUM_VEC already set to 1, so no else necessary
    
	
	
    //for( j=0; j<ktn; j++ ) { Rprintf( "%g \t" , BIG_SUM_VEC[j] ) ; }
    
	for( i=0; i<n; i++ ) { ////// for each site
		
        i_x_tau = i * tau ;
		
		if( geodesic[0] == 1 ) {
			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				n_vec_prod = A_cc[i]*B_cc[j] + A_cs[i]*B_cs[j] + A_s[i]*B_s[j] ;
				if(  n_vec_prod < -1 ) { n_vec_prod = -1 ; }
				if(  n_vec_prod > 1 ) { n_vec_prod = 1 ; }
				term1 = 6378.1 * acos( n_vec_prod ) ;
				term2 = pow(  pow( term1, 2 )  +   pow( A_z[i] - B_z[j], 2 ) , 0.5 ) ;
				w_vec[ j ] = term2 * BIG_SUM_VEC[j] ;
				//Rprintf( "%g \t" , w_vec[ j ] );
				//Rprintf( "%g \t" , term1 );
			}
			
		} else {
			
			for( j=0; j<ktn; j++ ) { //// for each site-lag
				term1 = pow(  pow( A_cc[i] - B_cc[j], 2 )  +  pow( A_cs[i] - B_cs[j], 2 )  +  pow( A_z[i] - B_z[j], 2 ),   0.5 ) ;
				term2 = term1 ;
				w_vec[ j ] = term2 * BIG_SUM_VEC[j] ;
				//Rprintf( "%g \t" , w_vec[ j ] );
				//Rprintf( "%g \t" , term1 );
			}
			
			
		}
		
		
        //Rprintf( "\n" );
		
        for( j=0; j<ktn; j++ ) {
            //w_vec[ j ] = pow( 2.71828183, malpha * w_vec[ j ] ) ;
			w_vec[ j ] = exp( malpha * w_vec[ j ] ) ;
            //Rprintf( "%g \t" , w_vec[ j ] );
        }
        //Rprintf( "\n" );
        
        if( start_self_ref[0] > -1 ) {
            for( j=0; j<len_ssr[0]; j++ ) { ///// start_self_ref must be zero-base vector
                w_vec[ start_self_ref[j]+i ] = 0 ;
                //Rprintf( "%g \t" , w_vec[ j ] );
            }
        }
        //Rprintf( "\n" );
        
        
        BIG_SUM_1 = 0 ;
        for( j=0; j<ktn; j++ ) {
            BIG_SUM_1 = BIG_SUM_1 + w_vec[ j ] ;
        }
        //Rprintf( "%g \t" , BIG_SUM_1 );
        
        if( BIG_SUM_1 == 0 ) { inv_wstnd = 1 ; } else { inv_wstnd = pow( BIG_SUM_1, -1 ) ; }
        //inv_wstnd = pow( BIG_SUM_1, -1 ) ;
		//Rprintf( "%g \t" , BIG_SUM_1 );
        
        inv_wstnd_x_wflatten = inv_wstnd * wflatten ;
        for( j=0; j<ktn; j++ ) {
            w_vec[ j ] = w_vec[ j ] * inv_wstnd_x_wflatten ;
            //Rprintf( "%g \t" , w_vec[ j ] );
        }
        //Rprintf( "\n" );
        
        snv_k = 0 ;
        for( j=0; j<ktn; j++ ) {
            if( !isnan(w_vec[ j ]) ) {
                if( w_vec[ j ] > sparse_cutoff ) {
                    sparse_ndx_vec[ snv_k ] = j ;
                    snv_k = snv_k + 1 ;
					//Rprintf( "%g \t" , w_vec[ j ] );
                }
				//Rprintf( "\n" );
            }
			
		}
        
        
        for( ii=0 ; ii<tau ; ii++ ) {
            BIG_SUM_2 = 0 ;
            for( j=0; j<snv_k; j++ ) {
				this_ndx = sparse_ndx_vec[j] ; 
				z_in_key = tau * z_rep_in[ this_ndx ] + ( z_lags_vec[ this_ndx ] + ii ) ;
				
				if( ( z_in_key >= 0 ) & ( z_in_key <= size_Zd ) ) {
					BIG_SUM_2 = BIG_SUM_2 + w_vec[ this_ndx ] * Z_in[ z_in_key ] ; /// can move flatten here
				}
            }
            Z_out[ ii + i_x_tau ] = BIG_SUM_2 ;
        }
		
	}
}





void distance_AB( double *A_x, double *A_y, double *B_x, double *B_y, double *D_out, int *d_A, int *d_B ) {
    int i, j ;
    int d_a = d_A[0] ;
    int d_b = d_B[0] ;
    for( i=0; i<d_a; i++ ) {
        for( j=0; j<d_b; j++ ) {
            D_out[ j*d_a + i ] = pow( pow( A_x[i] - B_x[j] , 2 ) + pow( A_y[i] - B_y[j] , 2 ), 0.5 ) ;
        }
    }
}




void distance_geodesic_AB( double *A_lat, double *A_lon, double *B_lat, double *B_lon, double *D_out, int *d_A, int *d_B ) {
    /// output MX will be d_a x d_b -- i.e., A is treated is X, B as Y in covariance Mx
    int i, j, d_a = d_A[0], d_b = d_B[0] ;
    double delta_lat , delta_lon, term1, term2 ;
    for( i=0; i<d_a; i++ ) {
        for( j=0; j<d_b; j++ ) {
            delta_lat = A_lat[i] - B_lat[j] ;
            delta_lon = A_lon[i] - B_lon[j] ;
            term1 = pow( sin(delta_lat*0.5), 2 ) ;
            term2 = pow( sin(delta_lon*0.5), 2 ) * cos(A_lat[i]) * cos(B_lat[j]) ;
            D_out[ j*d_a + i ] = 6378.1  *  2*asin( pow( term1+term2, 0.5 ) ) ;
        }
    }
}



//// cd /Users/dzes/Files/Creations/C ; R64 CMD SHLIB widals.c





//// cd /Users/dzes/Files/Creations/C ; R64 CMD SHLIB widals.c



