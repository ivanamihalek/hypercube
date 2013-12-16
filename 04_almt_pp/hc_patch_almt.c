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


typedef struct {
    int seqno;
    int aligned;
    double pct_id;
    double pct_sim;
} Dist_descr;


int patch_almt ( Alignment * alignment) {

    int i,j, k, no_seqs = alignment->number_of_seqs;
    int group_i, group_size, group_ctr, column_length;
    int chunk_size;
    int retval;
    char logname [BUFFLEN];
    int pos, patched, pos_patched;
    double pct_id, pct_sim;
    
    char *seq_i, *seq_j;
    int * int_arr_ptr;
    char ** new_sequence, **old_sequence;
    Dist_descr *nbr;
    FILE * log;

    //fprintf (stderr, "fix before using: (i) patching only within group (ii) no patching in qry\n");
    //exit (1);
    column_length = 	alignment->number_of_seqs;
    
    sprintf (logname, "%s.patchlog", options.outname);
    if ( ! (log = efopen( logname, "w"))  ) return 1;

    fprintf (log, "similarity cutoff: %8.3lf\n",  options.patch_sim_cutoff);
    fprintf (log, "length cutoff: %8.3lf\n",  options.patch_min_length);


    /*********************************/
    /*********************************/
    /* get rid of trivial positions: */
    int new_length = alignment->length;
    char * struct_seq = alignment->struct_seq;

    
    for (pos= alignment->length-1; pos >= 0; pos --) {
	if ( struct_seq  && struct_seq[pos] != '.' ) continue;
	if ( alignment->column_gaps[pos]+alignment->column_unks[pos] < column_length) continue;
	
	/* splice (clunker; fix at some other time)*/
	chunk_size = alignment->length - pos - 1;
	if ( chunk_size < 0 ) {
	    fprintf (stderr, "Error in %s:%d.\n", __FILE__, __LINE__);
	    exit (1);
	}
	i = 0;
	for (i=0; i<no_seqs; i++) {
	    seq_i = alignment->sequence[i];
	    memmove (seq_i+pos, seq_i+pos+1, chunk_size*sizeof(char));
	    alignment->seq_gaps[i] --;
	}
	seq_i = alignment->struct_seq;
	memmove (seq_i+pos, seq_i+pos+1, chunk_size*sizeof(char));
	int_arr_ptr = alignment->sunk;
	memmove (int_arr_ptr+pos, int_arr_ptr+pos+1, chunk_size*sizeof(int));

	int_arr_ptr = alignment->column_gaps;
	memmove (int_arr_ptr+pos, int_arr_ptr+pos+1, chunk_size*sizeof(int));
	int_arr_ptr = alignment->column_unks;
	memmove (int_arr_ptr+pos, int_arr_ptr+pos+1, chunk_size*sizeof(int));
	for (group_ctr = 0; group_ctr < alignment->no_groups; group_ctr++) {
	    int_arr_ptr = alignment->group[group_ctr].gaps;
	    memmove (int_arr_ptr+pos, int_arr_ptr+pos+1, chunk_size*sizeof(int));
	    int_arr_ptr = alignment->group[group_ctr].unks;
	    memmove (int_arr_ptr+pos, int_arr_ptr+pos+1, chunk_size*sizeof(int));
	}
	new_length--;
    }
    alignment->length = new_length;

    new_sequence = chmatrix (alignment->number_of_seqs, alignment->length);
    if ( !new_sequence ) return 1;

    /* the length has changed (perhaps) so we have to be a bit careful
       about copying to the new place */
    for (i=0; i<no_seqs; i++) {
	memcpy ( new_sequence[i], alignment->sequence[i], alignment->length*sizeof(char) );
    }
    old_sequence =  alignment->sequence;
    alignment->sequence = new_sequence;
    

    /************************************************************/
    /* get rid of strechces of sequence of suspicious content   */
    if (options.guess_bad_translation) {
	int find_bad_transl (Alignment *alignment, FILE * log);
	find_bad_transl (alignment, log);
    }


    /* seq similarity  */
    retval = seq_similarity_indicators (alignment);
    if (retval) return retval;

 
    /************************************************************/
    /* for each sequence, except query (the reference sequence), */
    /* sort the remaining sequences in the order of preference  */
    /* for patching */
    if ( ! (nbr=emalloc (no_seqs*sizeof (Dist_descr) )) ) return 1;

    for (i=0; i<no_seqs; i++) {
	
	group_i = alignment->belongs_to_group[i];
	/* never patch the representative  sequence ( the first sequence in a group)*/
	if ( i == alignment->group[group_i].group_member[0] ) continue;
	
	nbr[0].seqno = -1; nbr[0].pct_id = 0.0; nbr[0].pct_sim = 0.0;

	for (j=0; j<no_seqs; j++) {
	    if (i==j) continue;
	    if (group_i != alignment->belongs_to_group[j]) continue;
	    pct_id  = (double)alignment->identical_sites[j][i]/alignment->aligned_sites[i][j];
	    pct_sim = (double)alignment->similar_sites[j][i]/alignment->aligned_sites[i][j];

	    /* find place in the nbr array */
	    k = 0;
	    while ( k<no_seqs && pct_id <nbr[k].pct_id) k++;
	    while ( k<no_seqs && pct_sim<nbr[k].pct_sim) k++;

	    /* move no_seqs-1-k elements of the nbr array;
	       and write this sequence in position k */
	    if (  k < no_seqs) {
		if ( k < no_seqs-1 ) {
		    memmove ( nbr+k+1, nbr+k, (no_seqs-1-k)*sizeof(Dist_descr) ); 
		}
		nbr[k].seqno = j;
		nbr[k].pct_sim = pct_sim;
		nbr[k].pct_id  = pct_id;
		nbr[k].aligned = alignment->aligned_sites[i][j];
	    }
	}

	seq_i = old_sequence[i];
	fprintf (log, "%s\n", alignment->name[i]);
	patched = 0;
	group_size = alignment->group[group_i].no_members;

	for (pos=0; pos < alignment->length; pos ++) {
	  /* we are not patching postns which are not gappped, and not 'X' */
	    if (  seq_i[pos] != 'X' && seq_i[pos] != '.' && seq_i[pos] != '-' )
		continue;
	    /*    nor the positions which are gap for most of the members of the group */
	    if ( (double) alignment->group[group_i].gaps[pos]/group_size > 0.5 )
		continue;
	    /* patch from the closest relative possible */
	    pos_patched = 0;
	    for (k=0; k<no_seqs; k++) {
		 if (nbr[k].pct_sim < options.patch_sim_cutoff) break; /* the nbr's are sorted */
		 if (nbr[k].aligned < options.patch_min_length*alignment->length) continue;
		 j = nbr[k].seqno;
		 if ( j < 0) break; /* how and why was this supposed to happen? - no candidate */
		 seq_j = old_sequence[j];
		 if ( seq_j[pos] != 'X' && seq_j[pos] != '.' && seq_j[pos] != '-') {
		      new_sequence[i][pos] = seq_j[pos];
		      patched = 1;
		      pos_patched = 1;
		      fprintf (log, "\tpos %3d from %s (aligned length: %3d   sim: %6.3lf   id: %6.3lf)\n",
			       pos+1, alignment->name[j], nbr[k].aligned, nbr[k].pct_sim, nbr[k].pct_id);
		      break;
		 }
	    }
	    if (!pos_patched) new_sequence[i][pos] = 'X'; /* conservation score will ignore it */
	}
	if ( !patched )  fprintf (log, "\t not patched\n");

    }

    free_cmatrix ( old_sequence );
     
    /* recalc seq similarity  */
    retval = seq_similarity_indicators (alignment);
    if (retval) return retval;

    free (nbr);
    fclose (log);
    
    return 0;
}


# if 0
/* junkyard */
printf ("%2d (%s) -- %2d (%s) ",
	i, alignment->name[i], j,  alignment->name[j]);
printf ("a:%3d  s:%3d(%6.3lf)  i:%3d(%6.3lf)\n", alignment->aligned_sites[i][j],
	alignment->similar_sites[j][i],   nbr[k].pct_sim,
	alignment->identical_sites[j][i], nbr[k].pct_id);
# endif


/******************************************/
int find_bad_transl (Alignment *alignment, FILE * log) {
    /* #### !!!!! SOME HARDCODED STUFF  !!! ### */
    /* these things are not getting patched - I don't quite understande why */

    int i, g, s, pos, aa;
    int a_index;
    int bad_stretch_begin, bad_stretch_length;
    int recovery_length;
    int no_seqs = alignment->number_of_seqs;
    int no_groups    = alignment->no_groups;
    int no_positions = alignment->length;
    char *seq_i;
    double ***freq;
    Group *group;

    if ( ! (freq          = d3matrix(no_positions, no_groups, no_of_aa_types))) return 1;
    /****************************************************************/
    /****  counting the aa type frequencies                     *****/
    /****************************************************************/
    for (pos=0; pos< no_positions; pos++) {
   
	for (g=0; g<no_groups; g++) {
       
 	    group = alignment->group+g;
	    //memset (freq[pos][g], 0, no_of_aa_types*sizeof(double));
	    for (a_index=0; a_index<no_of_aa_types; a_index++) {
		 freq[pos][g][a_index] = 0.0;
	    }

	    for (i=0; i< group->no_members; i++) {
		s  = group->group_member[i];
		aa = (int) alignment->sequence[s][pos];
		/* skip x when counting frequencies  */
		if ( aa == 'X' || aa == '.' || aa == '-') continue;
		a_index = aa2index[aa];
       
		freq[pos][g][a_index] ++;
	    }
	    
	    normalize (freq[pos][g], no_of_aa_types, 1);

	}
    }
    
    for (i=0; i<no_seqs; i++) {
	
	seq_i = alignment->sequence[i];
	g = alignment->belongs_to_group[i];
	
	bad_stretch_length = 0;
	bad_stretch_begin = -1;
	recovery_length = 0;
	for (pos=0; pos< no_positions; pos++) {
	    
	    aa = (int)seq_i[pos];
	    if ( aa == 'X' || aa == '.' || aa == '-') continue;
	    
	    a_index = aa2index[aa];
	    if ( freq[pos][g][a_index] < 0.1) {
		
		if ( bad_stretch_length == 0 ) bad_stretch_begin = pos;
		bad_stretch_length ++;
		recovery_length = 0;
		
	    } else if ( bad_stretch_begin > 0 ) {
		
		recovery_length ++;
		if ( recovery_length > 2  ||  bad_stretch_length  < 5 ) {
		    
		    
		} else {
		    int pos_x;
		    for (pos_x = bad_stretch_begin; pos_x < pos; pos_x++) {
			seq_i[pos_x] = 'X';
		    }
		    fprintf (log, " %s has a bad stretch from %d to %d (will attempt patching)\n", 
			     alignment->name[i], bad_stretch_begin, pos-1);
		}
		
		bad_stretch_begin = -1;
		bad_stretch_length = 0;
	    }
	}
    }

    free_d3matrix(freq);


    return 1;
}
