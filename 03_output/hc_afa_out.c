# include "hypercube.h"

int afa_out (Options * options, Alignment * alignment) {

    int i, pos;
    char * seq;
    FILE * fptr;
    char filename[BUFFLEN];
 
    sprintf (filename, "%s.patched.afa", options->outname);

    fptr = efopen (filename, "w");
    if (!fptr) return 1;
    
    for (i=0; i<alignment->number_of_seqs; i++) {
	
	fprintf ( fptr,  ">%s\n", alignment->name[i]);
	seq = alignment->sequence[i];
	
	for (pos=0; pos<alignment->length; pos++) {
	    if ( seq [pos] == '.' ) {
		/* seaview doesn't like dots */
		fprintf ( fptr,  "-");
	    } else {
		fprintf ( fptr,  "%c",  seq [pos]);
	    }
	    if ( pos%50 == 49 ) fprintf (fptr, "\n");
	}
	if ( pos%50 ) fprintf (fptr, "\n");
  
    }
    fclose (fptr);
    
    return 0;
}
