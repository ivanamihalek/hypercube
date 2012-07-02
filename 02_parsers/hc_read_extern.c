# include "hypercube.h"

int read_extern_spec_method ( char *filename,  Alignment * alignment, double *spec_score) {

    FILE * fptr;
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    char comment_char, type, orig_type;
    int line_ctr, retval;
    int max_token;
    int pos, rep_seq;
   
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;

    printf ("reading %s.\n", filename);
    
    line_ctr = 0;
    memset ( line, 0, LONGSTRING);
    while(fgets(line,LONGSTRING,fptr)!=NULL){
 	line_ctr++;
 	/* tokenize */
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	switch ( retval ) {
	case  TOK_TOOMNY:
	    fprintf (stderr, "In %s: too many tokens on line %d.", filename, line_ctr);
	    exit (1);
	    break;
	case TOK_TOOLONG:
	    fprintf (stderr, "In %s: token too long on line %d.", filename, line_ctr);
	    exit (1);
	    break;
	}
	if ( max_token < 0 ) continue;
	if ( max_token < 2 ) {
	    fprintf (stderr, "In %s: too few tokens on line %d.", filename, line_ctr);
	    exit (1);
	}
	pos   = atoi (token[0]);
	if ( pos < 0 || pos > alignment->length) {
	    fprintf (stderr, "In %s: negative position or pos out of range: line %d.", filename, line_ctr);
	    exit (1);
	}
	type  = token[1][0];
	if (type == '-') type = '.';
	rep_seq   = alignment->group[0].group_member[0];
	orig_type = alignment->sequence[rep_seq][pos];
	if ( orig_type != type ) {
	    fprintf (stderr, "Type mismatch in %s: %c (vs %c in the original alignment).\n",
		     filename, type, orig_type);
	}
	spec_score[pos] = atof (token[2]);
    }


    fclose (fptr);

    
 

    
    
    return 0;
}
