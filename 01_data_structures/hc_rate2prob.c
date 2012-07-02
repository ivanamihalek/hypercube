# include "hypercube.h"
// Felsenstein, chap 13, eq 13.15 through 13.19

# define N 20


int rate2prob (double ** rate_sym, double * freq, double time_step, int no_timesteps, double *** prob_matrix) {

    int i, j, k, clock;
    double sum = 0;
    double VL[N][N], VR[N][N], egvl[N];
    int factorize (double ** rate_sym, double * freq, double VL[N][N],double VR[N][N],double egvl[N]);

    /* factorization */
    factorize (rate_sym, freq, VL, VR, egvl);
    
   
   
    /* the exponent: */
    for (clock=0; clock < no_timesteps; clock++ ) {
	 for (i=0; i<N; i++) {
	     for (j=0; j<N; j++) {
		 prob_matrix[clock][i][j] = 0;
		 for (k=0; k<N; k++) {
		     prob_matrix[clock][i][j] += VL[i][k]*exp(egvl[k]*clock*time_step)*VR[k][j];
		 }
	     }
	 }
	 /* sanity check: is this a prob matrix? */
	 for (j=0; j<N; j++) {
	     sum = 0.0;
	     for (i=0; i<N; i++) {
		 sum +=   prob_matrix[clock][i][j];
	     }
	     if ( fabs (1.0-sum) > 1.e-10) {
		 fprintf (stderr, "Numerical in rate2prob()? \n");
		 fprintf (stderr, "clock %3d   column %2d:  fabs(1.0-column_sum)=%8.1le\n",
			  clock, j, fabs (1.0-sum));
		 exit (1);
	     }
	 }
    }

    

# if 0
    /* at very long times,  each column should be equal to stationary freq */
     clock = no_timesteps-1;
     //clock = (int)(clock*0.75);
     printf ("\n");
 
     for (i=0; i<N; i++) {
	 for (j=0; j<N; j++) {
	     printf ("%6.3lf ", prob_matrix[clock][i][j] );
	 }
	 printf (" **%6.3lf \n", freq[i]);
     }
     printf ("\n");

     exit (1);
# endif
    
     return 0;

}



/******************************************************************************************/

int factorize (double ** rate_sym, double * freq, double VL[N][N],double VR[N][N],double egvl[N]) {

    int i, j, k;
    double sum = 0, norm;
    double rate[N][N];
    double a[N][N];
    double b[N][N];
    double c[N][N];
    double d[N][N];
 
    double **A;
    int n = N, lda = N, info, lwork;
    double wkopt;
    double* work;
    void dsyev_ ( char* jobz, char* uplo, int* n, double* a, int* lda,
		  double* egvl, double* work, int* lwork, int* info);
    
    
    if (! ( A=dmatrix(N,N))) return 1;
 
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    if (i==j) continue;
	    rate[i][j] = freq[i]*rate_sym[i][j];
	}
    }

    /* normalize  */
    norm = 0;
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    if (i==j) continue;
	    norm += rate[i][j]*freq[j];
	}
    }

    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    if (i==j) continue;
	    rate[i][j]     /= norm;
	    rate_sym[i][j] /= norm;
	}
    }


    /* the columns should add up to zero */
    for (i=0; i<N; i++) {
	sum = 0.0;
	for (j=0; j<N; j++) {
	    if (i==j) continue;
	    sum += rate[j][i];
	}
	rate [i][i] = -sum;
    }
   

    /* write rate as a product of a diagonal and sym matrix */
    for (i=0; i<N; i++) {
	 
	for (j=i+1; j<N; j++) {
	    b[i][j] = b[j][i] = rate_sym[i][j];
	    d[i][j] = d[j][i] = 0.0;
	}
	 
	b[i][i]  = 0;
	for (j=0; j<N; j++) {
	    if (i==j) continue;
	    b[i][i] -= freq[j]* b[i][j];
	}
	b[i][i] /= freq[i];
	d[i][i]  = freq[i];
    }

    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    a[i][j] = 0;
	    for (k=0; k<N; k++) {
		a[i][j] += d[i][k]*b[k][j];
	    }
	}
    }

    /* c, our new symmetric matrix to be diagonalized */
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    c[i][j] = c[j][i] = b[i][j]*sqrt(freq[i]*freq[j]);
	}
    }
  

    /* transpose the matrix c  or just, copy, doesn't matter*/
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    A[i][j] = c[j][i];
	}
    }

    lwork = -1;
    dsyev_ ( "Vectors", "Upper", &n, A[0], &lda, egvl, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work  = (double*)malloc( lwork*sizeof(double) );
    /* Solve eigenproblem */
    dsyev_ ( "Vectors", "Upper", &n, A[0], &lda, egvl, work, &lwork, &info );

    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    VL[i][j] = A[j][i]*sqrt(freq[i]);
	    VR[i][j] = A[i][j]/sqrt(freq[j]);
	}
    }


    
 # if 0 
     for (i=0; i<N; i++) {
	 for (j=0; j<N; j++) {
	     a[i][j] = 0;
	     for (k=0; k<N; k++) {
		 a[i][j] += VL[i][k]*exp(egvl[k])*VR[k][j];
	     }
	 }
     }
     
     for (i=0; i<N; i++) {
	 for (j=0; j<N; j++) {
	     printf ("  %6.3lf", a[i][j]);
	 }
	 printf ("\n");
     }
     printf ("check:\n");
     for (i=0; i<N; i++) {
	 sum = 0.0;
	 for (j=0; j<N; j++) {
	     sum +=  a[j][i];
	 }
	 printf (" %3d  %6.3lf\n", i , sum);
     }
   
     exit (1);
# endif

     free_dmatrix(A);
     free(work);
     return 0;

}
