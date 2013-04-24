# include "hypercube.h"

int no_kwds = NO_KWDS; /* unless it turns out to be smaller, in set_hwds() below */

int  read_opt_file (char * filename ) {
    
    FILE * fptr, *log = stdout;
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    char comment_char;
    int ctr, line_ctr, retval, token_ctr;
    int max_token, token_assigned;
    /***************/
    int echo_kwds ();
    int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
		 char * fmt, char * warnstr) ;
    /***************/
    
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;

    
    line_ctr = 0;
    memset ( line, 0, LONGSTRING);
    while(fgets(line,LONGSTRING,fptr)!=NULL){
 	line_ctr++;
 	/* tokenize */
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	switch ( retval ) {
	case  TOK_TOOMNY:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Too many tokens.");
	    fclose (log);
	    break;
	case TOK_TOOLONG:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Token too long.");
	    fclose (log);
	    break;
	}
	if ( max_token < 0 ) continue;
	
	
	token_assigned = 0;
	/* check first token for meaning (first token should be a keyword)*/
	for (ctr=0; ctr<no_kwds && ! token_assigned; ctr++) {
	    
	    if (  strcmp (token[0], kwd[ctr].name)) continue;
	    
	    if ( kwd[ctr].type != OPTN_FLAG  && max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\tKeyord %s should be followed by value.\n",
			 token[0]);
		return 1;
	    }

		
	    switch (kwd[ctr].type) {
		
	    case OPTN_FLAG:
		/*toggle*/
		*(int *)kwd[ctr].storage = 1 - *(int *)kwd[ctr].storage;
		break;
		
	    case OPTN_INTEGER:
		*(int *)kwd[ctr].storage = atoi( token[1] );
		break;
		
	    case OPTN_DOUBLE:
	    case OPTN_EXP:
		*(double *)kwd[ctr].storage = atof( token[1] );
		break;
		
		
	    case OPTN_STRING:
		sprintf ( kwd[ctr].storage, "%s", token[1]);
		break;	

	    case OPTN_INT_ARRAY:
		for (token_ctr=1; token_ctr <= max_token; token_ctr++) {
		    *((int *)kwd[ctr].storage+token_ctr)  = atoi (token[token_ctr]);
		}
		*(int*)kwd[ctr].storage = max_token; /*store the array size here*/
		break;
		
	    case OPTN_STRING_ARRAY:
		if (max_token >= MAX_STRING_ARRAY) {
		    fprintf (stderr,
			     "Error reading %s: too many tokens. (Increase MAX_CHAR_ARRAY and recompile).\n",
			     token[0]);
		    return 1;
		}
		if (max_token < 2 || max_token%2 ) {
		    errmsg ( log, line_ctr, line,
			 "\tKeyord %s should be followed by filename and tag (repeated for each extern method).\n",
			     token[0]);
		    return 1;
		}
		char ***string_array = (char ***)kwd[ctr].storage; /* gross, but can't help it right now */
		*string_array =  chmatrix (MAX_STRING_ARRAY, BUFFLEN);
		if ( ! *string_array ) {
		    fprintf (stderr, "Error allocating space for %s.\n", token[0]);
		    return 1;
		} 
		for (token_ctr=1; token_ctr <= max_token; token_ctr++) {
		    sprintf ( (*string_array)[token_ctr], "%s", token[token_ctr]);
		}
		(*string_array)[0][0] = max_token; /*store the array size here*/
		break;
		
		
	    default:
		fprintf ( stderr, "Unrecognized kwd type: %s.\n", token[0]);
		return 1;
	    }
	    
	    token_assigned = 1;
	}
	
	if ( ! token_assigned){
	    errmsg  (stderr, line_ctr, line, "\t Keyword %s not recognized.\n", token[0]);
	    fprintf (stderr, "Recognized keywords (followed by default values):\n");
	    echo_kwds();
	    return 1;
	}    
	
	memset (line, 0, LONGSTRING);
   }
   fclose (fptr);

   /* check fot stts option: */
   if ( options.pdbname[0] && ! options.struct_name[0]) {
       fprintf (stderr,"%s%s%s%s", 
		"pdb file (", options.pdbname, ") must be accompanied by the corresponding sequence ",
		"in the alignment (\"pdb_name\" kwd)\n");
       return 1;
   }
   
   /* another bit of compatibility */
   if ( options.compare && !options.rate_matrix_file[0]) {
       fprintf (stderr,"kwd \"compare\" must be accompanied by the name of the rate matrix file.\n");
       return 1;
   }
   
   echo_kwds ();
   
   return 0;
}

/****************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
	     char * fmt, char * warnstr) {

    fprintf ( fptr, "Error on line %3d:     %s", line_ctr, line);
    fprintf ( fptr, fmt, warnstr);
    return 0;
}
/****************************************************************************/
int echo_kwds () {
    int ctr;

    for (ctr=0; ctr<no_kwds; ctr++) {
	switch (kwd[ctr].type) {
		
	case OPTN_FLAG:
	    if (  *(int *)kwd[ctr].storage ) {
		fprintf ( stderr, "\t%30s     ON\n",
			  kwd[ctr].name);
	    } else {
		fprintf ( stderr, "\t%30s     OFF\n",
			  kwd[ctr].name);
	    }
	    break;
		
	case OPTN_INTEGER:
	    fprintf (stderr, "\t%30s %8d\n", kwd[ctr].name, *(int *)kwd[ctr].storage);
	    break;
		
	case OPTN_DOUBLE:
	    fprintf (stderr, "\t%30s %8.2lf\n", kwd[ctr].name, *(double *)kwd[ctr].storage);
	    break;
		
	case OPTN_EXP:
	    fprintf (stderr, "\t%30s %12.2le\n", kwd[ctr].name, *(double *)kwd[ctr].storage);
	    break;
		
	case OPTN_STRING:
	    if (  *(char *)kwd[ctr].storage ) {
		fprintf ( stderr, "\t%30s     %s\n",
			  kwd[ctr].name, (char *)kwd[ctr].storage);
	    } else {
		fprintf (stderr, "\t%30s     (none)\n", kwd[ctr].name);
	    }
	    break;

	case OPTN_INT_ARRAY:
	    fprintf (stderr, "\t%30s     %s\n", kwd[ctr].name, "(int array)");
	    break;
	    
	case OPTN_STRING_ARRAY:
	    fprintf (stderr, "\t%30s     %s\n", kwd[ctr].name, "(string array)");
	    break;
	    
		
	default:
	    fprintf (stderr, "Unrecognized kwd type (oink!?)\n");
	    return 1;
	}
    }
    return 0;

}
/************************************/
/************************************/
/* this function is set in main, so the keywords can
   be echoed by "options" on  the commandline */
int set_keywords () {

    int ctr = -1;
    memset (kwd, 0, NO_KWDS*sizeof(Keyword));
    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "acc");
    kwd[ctr].type    = OPTN_DOUBLE;
    kwd[ctr].storage = &(options.acc_cutoff);
    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "almtname");
    kwd[ctr].type = OPTN_STRING;
    kwd[ctr].storage = options.almtname;

    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "compare");
    kwd[ctr].type = OPTN_FLAG;
    kwd[ctr].storage = &(options.compare);
    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "dssp");
    kwd[ctr].type    = OPTN_STRING;
    kwd[ctr].storage = options.dssp_name;
    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "exchangeability");
    kwd[ctr].type = OPTN_FLAG;
    kwd[ctr].storage =  &(options.exchangeability);
    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "extern");
    kwd[ctr].type    = OPTN_STRING_ARRAY;
    kwd[ctr].storage = &(options.extern_spec_methods);
    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "groups");
    kwd[ctr].type    = OPTN_STRING;
    kwd[ctr].storage = options.groups;

    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "max_gaps");
    kwd[ctr].type    = OPTN_DOUBLE;
    kwd[ctr].storage = &(options.max_gaps);

    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "outname");
    kwd[ctr].type    = OPTN_STRING;
    kwd[ctr].storage = options.outname;
 
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "patch_sim_cutoff");
    kwd[ctr].type    = OPTN_DOUBLE;
    kwd[ctr].storage = &(options.patch_sim_cutoff);

    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "patch_min_length");
    kwd[ctr].type    = OPTN_DOUBLE;
    kwd[ctr].storage = &(options.patch_min_length);

    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "guess_bad_transl");
    kwd[ctr].type    = OPTN_FLAG;
    kwd[ctr].storage =  &(options.guess_bad_translation);
    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "pdb_file");
    kwd[ctr].type = OPTN_STRING;
    kwd[ctr].storage = options.pdbname;

    
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "raw");
    kwd[ctr].type = OPTN_FLAG;
    kwd[ctr].storage = &(options.raw_score);
        
 
    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "rate_matrix");
    kwd[ctr].type = OPTN_STRING;
    kwd[ctr].storage = options.rate_matrix_file;


    /*******************/
    if ( (++ctr) >= NO_KWDS ) return 1;
    sprintf (kwd[ctr].name, "%s",  "pdb_name");
    kwd[ctr].type = OPTN_STRING;
    kwd[ctr].storage = options.struct_name;

    /* this guy is global: */
    no_kwds = ctr+1;
 
 
    
    return 0;
}
