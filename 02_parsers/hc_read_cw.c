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

int read_clustalw (Alignment * alignment){
    
    FILE * fptr = NULL;
    char line[BUFFLEN];
    int  number_of_seqs, almt_length, ctr, ctr2;
    int * seq_pos, pos,struct_pos_ctr;
    int  struct_found;
    int * pos_ctr;
    char * seq_ptr;
    char ** sequence;
    char ** name;
    char curr_name[BUFFLEN];
     
    /* open file */
    fptr = efopen ( options.almtname, "r");
    if ( !fptr ) return 1;

    memset (alignment, 0, sizeof(Alignment) );
    
     
    /* find the alignment length info */
    almt_length = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "MSF:" ) ){
	    sscanf (line, "%*s %d", &almt_length);
	    break;
	}
    }
    if ( almt_length ) {
	/* printf ( "Alignment length in %s is %d.\n", cwname, almt_length); */
    } else {
	fprintf ( stderr, "Alignment length info not found in %s. Is the format gcg?\n",
		  options.almtname );
	return 1;
    }

    /* determine the number of sequences */
    number_of_seqs = 0;
    struct_found = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( ! strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    sscanf (line, "%*s %s", curr_name);
	    if (!strcmp(curr_name, options.struct_name)) {
		struct_found = 1;
	    } else {
		number_of_seqs++;
	    }
	}
    }
    if ( number_of_seqs ) {
	printf ( "Number of sequences in %s is %d.\n",  options.almtname, number_of_seqs);
    } else {
	fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n",
		  options.almtname);
	return 1;
    }

    if ( options.struct_name[0] && !struct_found ) {
 	fprintf ( stderr, "Structure  %s not  found in %s.\n",
		  options.struct_name,  options.almtname);
	return 1;
   }
    
    
    /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;
        
    alignment->struct_seq = emalloc (almt_length*sizeof(char));
    if ( !alignment->struct_seq ) return 1;
    
    alignment->sunk = emalloc (almt_length*sizeof(int));
    if ( !alignment->sunk ) return 1;
    
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;

    seq_pos = (int *) emalloc ( number_of_seqs*sizeof(int));
    if ( !seq_pos ) return 1;

    /****************************/
    /* read in                  */
    rewind(fptr);
    ctr = 0;
    struct_pos_ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if (!  strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    sscanf (line, "%*s %s", curr_name);
	    if ( !strcmp (curr_name, options.struct_name)) continue;
	    sprintf (name[ctr], "%s", curr_name);
	    ctr ++;
	}
    }
    
    /* check for duplicate names */
    for (ctr = 0; ctr <number_of_seqs;  ctr++) {
	for (ctr2 = ctr+1; ctr2 <number_of_seqs;  ctr2++ ) {
	    if ( ! strcmp (name[ctr], name[ctr2]) ) {
		fprintf ( stderr, "Duplicate names  found in the header of %s:  %s (names %d and %d)\n",
			  options.almtname, name[ctr], ctr, ctr2);
		return 1;
	    }
	}
    }
    
   
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( isspace (line[0] ) ) continue;
	sscanf (line, "%s", curr_name);

	
	if ( !strcmp (curr_name, options.struct_name)) {
	    seq_ptr = alignment->struct_seq;
	    pos_ctr = &struct_pos_ctr;
	    
	    
 	} else {
	    ctr = 0;
	    while (  ctr <number_of_seqs &&  strcmp (name[ctr], curr_name) ) ctr++;
	    if ( ctr >= number_of_seqs ) {
		fprintf ( stderr, "The name %s not found in the header of %s.\n",
			  curr_name,  options.almtname);
		return 1;
	    }
	    seq_ptr = sequence [ctr];
	    pos_ctr = seq_pos + ctr;
	}
	pos = 0;
	while ( ! isspace(line[pos]) ) pos++;
	while  (line[pos] != '\n' && pos < BUFFLEN) {
	    if ( !  isspace(line[pos] ) ){
                /* --> turn to uppercase */
		if ((line[pos]>=97) && (line[pos]<=122)) {line[pos] -= 32;}
		/* turn dash to dot */
		if ( line[pos]==45 )                   {line[pos]  = 46;} 
		/* turn tweedle to dot */
		if ( line[pos]==126)                   {line[pos]  = 46;} 
		/* turn X to dot */
		//if ( line[pos]==88)                    {line[pos]  = 46;}
		
		seq_ptr [ *pos_ctr ] = line[pos];
		(*pos_ctr)++;
	    }
	    pos ++;
	}
    }
    fclose(fptr);

    /* sanity check */
    for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	if ( seq_pos[ctr] >  almt_length ) {
	    fprintf (stderr,
		     "Sequence %s is longer (%d positions) than the alignment (%d positions).\n",
		     name[ctr],  seq_pos[ctr], almt_length);
	    return 1;
	} else if ( seq_pos[ctr] <  almt_length ) {
	    fprintf (stderr,
		     "Sequence %s is shorter (%d positions) than the alignment (%d positions).\n",
		     name[ctr],  seq_pos[ctr], almt_length);
	    return 1;
	}
    }

     
    alignment->number_of_seqs = number_of_seqs;
    alignment->length         = almt_length;
    alignment->sequence       = sequence;
    alignment->name           = name;

    /* free */
    free (seq_pos);

    
    return 0;
}
