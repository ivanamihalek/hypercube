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


int read_rate_matrix (char * infile_name, double **rate_sym, double *freq) {
    
    char comment_char;
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    int  max_token, tok_ctr, retval;
    int line_ctr;
    int i, j, ctr;
    FILE * fptr;
    int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING], char * fmt, char * warnstr) ;

    /* open table file */
    if ( !(fptr = efopen (infile_name, "r" )) ) return 1;
    memset ( line, 0, LONGSTRING);
    line_ctr = 0;
    ctr = 1;
    while(fgets(line,LONGSTRING,fptr)!=NULL){
	line_ctr++;
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	switch ( retval ) {
	case  TOK_TOOMNY:
	    errmsg ( stderr, line_ctr, line, "\t\t %s\n", "Too many tokens.");
	    return 1;
	    break;
	case TOK_TOOLONG:
	    errmsg ( stderr, line_ctr, line, "\t\t %s\n", "Token too long.");
	    return 1;
	    break;
	}
	if ( max_token < 0 ) continue;
	
	if ( ctr >  20) {
	    fprintf ( stderr, "Error reading blosum tables (line %d) - too many lines.\n", line_ctr);
	    return 1;
	}
	if ( max_token >= 20) {
	    fprintf ( stderr, "Error reading blosum tables (line %d) - too many tokens.\n", line_ctr);
	    return 1;
	}
	if ( ctr ==  20) {
	    for (tok_ctr=0; tok_ctr <= max_token; tok_ctr ++ ) {
		freq[tok_ctr] = atof(token[tok_ctr]);
	    }
	} else {
	    for (tok_ctr=0; tok_ctr <= max_token; tok_ctr ++ ) {
		rate_sym[ctr][tok_ctr] = atof(token[tok_ctr]);
	    }
	}
	ctr++;
   }
    
    fclose (fptr);
    
    for(i=0;i<20;i++){
	for (j=i+1;j<20;j++){
	    rate_sym[i][j] = rate_sym[j][i];
	}
    }

# if 0
    for(i=0;i<20;i++){
	for (j=0;j<20;j++){
	    printf ("%8.3lf", rate_sym[i][j]);
	}
	printf ("\n");
    }
    printf ("\n");
    for (j=0;j<20;j++){
	printf ("%8.3lf", freq[j]);
    }
    printf ("\n");
    exit (1);
# endif
    
    return 0;
}

/*********************************************************************************************/
