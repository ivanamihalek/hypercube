/*
This source code is part of specs, an application for
protein residue specialization and conservation scoring.
Written by Ivana Mihalek.
Copyright (C) 2010-2013 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# include "hypercube.h"


/***********************************************************************/
int chisquare ( int * population_1, int * population_2, int no_bins, int kstr, double *df, double *chsq, double *prob) {

    /* normalize distributions */

    int gammq(double a, double x, double * retval);
    int j;
    double total_1, total_2;
    double temp;
    double p_prefactor, n_prefactor;
    
    total_1 = 0;
    total_2 = 0;
    for (j=0; j < no_bins;j++) {
	total_1 +=  population_1[j];
	total_2 +=  population_2[j];
    }

    p_prefactor = sqrt (total_2/total_1);
    n_prefactor = 1.0/p_prefactor;
  
    *df = no_bins - kstr;
    *chsq = 0.0;
    for (j=0; j < no_bins;j++) {
	if (population_1[j] == 0.0 && population_2[j] == 0.0) {
	    --(*df);   
	} else {  
	    temp   = p_prefactor* population_1[j] - n_prefactor*population_2[j];
	    *chsq += temp*temp/(population_1[j]+population_2[j]);
	}
    }
    if ( gammq(0.5*(*df),0.5*(*chsq), prob) ) return 1;
    

    //printf ( " %8.1le  %8.1le     %8.1le \n", *df, *chsq, *prob);


    return 0;
}
/*************************************************************/
//Returns the incomplete gamma function P (a, x).
double gammp(double a, double x)
{
    int gcf(double *gammcf, double a, double x, double *gln);
    int gser(double *gamser, double a, double x, double *gln);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) {
	fprintf (stderr, "Invalid arguments in routine gammp (x=%5.1lf, a=%5.1lf).\n", x, a);
	exit (1);
    }
    if (x < (a+1.0)) {   // Use the series representation.
	gser(&gamser,a,x,&gln);
	return gamser;
    } else {   // Use the continued fraction representation
	gcf(&gammcf,a,x,&gln);
	return 1.0-gammcf;    //and take its complement.
    }
}
/*************************************************************/
//Returns the incomplete gamma function Q(a, x)  1 - P (a, x).
int gammq(double a, double x, double * retval) {
    
    int gcf(double *gammcf, double a, double x, double *gln);
    int gser(double *gamser, double a, double x, double *gln);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) {
	fprintf (stderr, "Invalid arguments in routine gammp (x=%5.1lf, a=%5.1lf).\n", x, a);
	exit (1);
    }
    if (x < (a+1.0)) {   // Use the series representation
	if (gser(&gamser,a,x,&gln)) {
	    return 1;
	}
	*retval = 1.0-gamser;   // and take its complement.
    } else {    //Use the continued fraction representation.
	if ( gcf(&gammcf,a,x,&gln) ) {
	    return 1;
	}
	*retval = gammcf;
    }

    return 0;
}
/*************************************************************/
//Returns the incomplete gamma function P (a, x) evaluated by its
//series representation as gamser.
//Also returns ln (a) as gln.
#define ITMAX 10000   //Maximum allowed number of iterations.
#define EPS 3.0e-7   //Relative accuracy.
int gser(double *gamser, double a, double x, double *gln)
{
    double gammln(double xx);
    int n;
    double sum,del,ap;
    *gln=gammln(a);
    if (x <= 0.0) {
	if (x < 0.0) {
	    fprintf ( stderr, "x less than 0 in routine gser\n");
	    return 1;
	}
	*gamser=0.0;
	return 0;
    } else {
	ap=a;
	del=sum=1.0/a;
	for (n=1;n<=ITMAX;n++) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (fabs(del) < fabs(sum)*EPS) {
		*gamser=sum*exp(-x+a*log(x)-(*gln));
		return 0;
	    }
	}
	fprintf( stderr, "a too large, ITMAX too small in routine gser\n");
	return 1;
    }
}
/*************************************************************/
//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction 
//representation as gammcf. Also returns ln (a) as gln.
#define FPMIN 1.0e-30   // Number near the smallest representable doubleing-point number.
int gcf(double *gammcf, double a, double x, double *gln)
{
    double gammln(double xx);
    int i;
    double an,b,c,d,del,h;
    *gln=gammln(a);
    b=x+1.0-a;    //Set up for evaluating continued fraction
    c=1.0/FPMIN;   // by modified Lentz's method (§5.2)
    d=1.0/b;    //with b0 = 0.
    h=d;
    for (i=1;i<=ITMAX;i++) {    //Iterate to convergence.
	an = -i*(i-a);
	b += 2.0;
	d=an*d+b;
	if (fabs(d) < FPMIN) d=FPMIN;
	c=b+an/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPS) break;
    }
    if (i > ITMAX) {
	fprintf(stderr, "a too large, ITMAX too small in gcf\n");
	return 1;
    }
    *gammcf=exp(-x+a*log(x)-(*gln))*h;   // Put factors in front.

    return 0;
}
/*************************************************************/
//Returns the value ln[(xx)] for xx > 0.
//Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
//accuracy is good enough.
double gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
/*************************************************************/
//Returns the complementary error function erfc(x) with fractional error everywhere less than
//1.2 × 10-7.
double erfcc(double x)
{
    double t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
			    t*(-0.82215223+t*0.17087277)))))))));
    // return (x >= 0.0) ? ans : 2.0-ans;
    return ans;
}
//Returns the complementary error function erfc(x).
double erffc(double x)
{
    double retval;
   
    if ( x < 0.0 ) {
	retval = gammp(0.5,x*x);
	return 1+ retval;
    } else {
	gammq(0.5,x*x, & retval);
	return retval;
    }
}

