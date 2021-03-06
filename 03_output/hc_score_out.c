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

/* this function assumes that spec_score contains ate least
   two scores; the first is assumed to be discriminant for all groups,
   and the second one determinant for the first group in groups file */

int  score_output (Alignment *alignment, Protein *protein, int *almt2prot, 
		   double ** score, double *** in_group_score, double ** spec_score){

    int almt_pos, score_ctr;
    int  pos, unsunk_pos, g, s;

    char filename[BUFFLEN];
    char aa = '\0', pdbid[PDB_ATOM_RES_NO_LEN+1] = {'\0'};
    double cvg;
    double ** cons_cvg = NULL;
    double ** spec_cvg = NULL;
    FILE  * fptr = NULL;
    Group * group;
    int nongap_pos_ctr[alignment->no_groups];
    int nongap_not_all_pos_ctr[alignment->no_groups];
   
    
    sprintf (filename, "%s.score", options.outname);

    fptr = efopen (filename, "w");
    if ( !fptr) return 1;

    /* find residue ranks (coverage)*/
    if ( !(cons_cvg = dmatrix(1+alignment->no_groups, alignment->length) )) return 1;
    if ( !(spec_cvg = dmatrix(1+alignment->no_groups, alignment->length) )) return 1;
	
    /* entropy */
    coverage (alignment, score[1], cons_cvg[0], -1);
    /* rvet */
    /* coverage (alignment, score[0], cons_cvg[0], -1); */
    
    /* conservation in the individual groups */
    for  (g=0; g < alignment->no_groups; g++) {
	coverage (alignment, in_group_score[0][g], cons_cvg[g+1], g+1); /* conservation */
	coverage (alignment, spec_score[g+1], spec_cvg[g+1], -(100+g)); /* discriminants */
    }
    /* spec score */
    g = 0;
    coverage (alignment, spec_score[g], spec_cvg[g], 0); /* determinants */
    


    /******** header **************************/
    fprintf (fptr, "%%%6s %8s", "almt", "gaps(%)");

    fprintf (fptr, " %8s ", "entropy");
    
    for  (g=0; g < alignment->no_groups; g++) {
	fprintf (fptr, " %s ", alignment->group[g].group_name);
    }
    
    fprintf (fptr, " %8s ", "discr");
    for (g=0; g<alignment->no_groups; g++) {
	fprintf (fptr, " dets_in_%s ", alignment->group[g].group_name);
    }
    
    for (g=0; g<alignment->no_groups; g++) {
	group = alignment->group+g;
	s  = group->group_member[0];
	fprintf ( fptr, " %s ", alignment->name[s]);
	fprintf ( fptr, " pos_in_%s ", alignment->name[s]);
    }
    
    
    if ( almt2prot ) {
	fprintf (fptr, " %8s ", "pdb_aa");
	fprintf (fptr, " %8s ", "pdb_id");
    }
    
    if ( options.dssp_name[0] ) {
      fprintf (fptr, "%6s", "surf");
    }

    fprintf (fptr, "\n");
   
    for (g=0; g<alignment->no_groups; g++) {
	nongap_pos_ctr[g] = -1;
	nongap_not_all_pos_ctr[g]  = -1;
    }
    pos = -1;
    unsunk_pos = -1;
    for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {

	if ( spec_score[0][almt_pos] != ALL_GAPS )  pos++;
	
	for  (g=0; g < alignment->no_groups; g++) {
	    group = alignment->group+g;
	    s  = group->group_member[0];
	    if (alignment->sequence[s][almt_pos] != '.' ) {
		nongap_pos_ctr[g]++;
		if ( spec_score[0][almt_pos] != ALL_GAPS) {
		    nongap_not_all_pos_ctr[g]++;
		}
	    }
	}

	fprintf (fptr, "%6d %8d", almt_pos+1, 
		 (int)(100*(double)alignment->column_gaps[almt_pos]/alignment->number_of_seqs));

	/************************/
	/* conservation scores  */
	/* overall conservation */
	score_ctr = 0; 
	if (alignment->sunk[almt_pos]) {
	    cvg = 1;
	} else {
	    unsunk_pos ++;
	    cvg = cons_cvg[score_ctr][unsunk_pos];
	}
	fprintf (fptr, " %8.2lf ", cvg);
	
	
	/*****************************/
	/* within-group conservation */
	for ( score_ctr=0; score_ctr<alignment->no_groups; score_ctr++) {
	    
	    group = alignment->group+score_ctr;
	    s  = group->group_member[0];
	    if (alignment->sequence[s][almt_pos] == '.' ) {//<========
		cvg = 1;
	    } else {
		cvg = cons_cvg[score_ctr+1][nongap_pos_ctr[score_ctr]];
	    }
	    fprintf (fptr, " %8.2lf ", cvg);
	    
	}
	/***********************/
	/* specificity  scores */
	score_ctr = 0;
	if ( spec_score[0][almt_pos] == ALL_GAPS  ||  alignment->sunk[almt_pos] ) {
	    cvg = 0.5;
	} else {
	    cvg = spec_cvg[score_ctr][ nongap_not_all_pos_ctr[score_ctr] ];
	}
	fprintf (fptr, " %8.2lf ", cvg);

	for (score_ctr = 1; score_ctr <= alignment->no_groups; score_ctr++) {

	    group = alignment->group+score_ctr-1;
	    s  = group->group_member[0];
	    if ( spec_score[0][almt_pos] == ALL_GAPS  ||  alignment->sequence[s][almt_pos] == '.' ) {
		cvg = 0.5;
	    } else {
		cvg = spec_cvg[score_ctr][ nongap_not_all_pos_ctr[score_ctr-1] ];
	    }
	    fprintf (fptr, " %8.2lf ", cvg);
	}

	
	
	/************************/
	/*  representative seqs */
	for (g=0; g<alignment->no_groups; g++) {
	    group = alignment->group+g;
	    s  = group->group_member[0];
	    aa = alignment->sequence[s][almt_pos];
	    fprintf (fptr, "%6c", aa);
	    if ( aa == '.' ) {
		fprintf (fptr, "%5s", ".");
	    } else {
		fprintf (fptr, "%5d", nongap_pos_ctr[g]+1);
	    }
	}
	if ( almt2prot) {
	     if ( almt2prot[almt_pos] >= 0 ) {
		  sprintf (pdbid, "%s", protein->sequence[ almt2prot[almt_pos] ].pdb_id );
		  aa = protein->sequence[ almt2prot[almt_pos] ].res_type_short;
		  fprintf (fptr, "%6c%5s", aa, pdbid);
	     } else {
		  fprintf (fptr, "%6s%5s", ".", ".");
	     }
	}
	
	/* surface accessibility */
	if ( options.dssp_name[0] ) {
	  int acc =  ( almt2prot[almt_pos] < 0 ) ?  
	    -1: protein->sequence[almt2prot[almt_pos]].solvent_accessible;
	    fprintf ( fptr," %6d",  acc);
	}
	
	fprintf (fptr, "\n");
    }

    fclose (fptr);
    free_dmatrix (cons_cvg);
    free_dmatrix (spec_cvg);
    
    return 0;

}


/************************************************************************/
int coverage ( Alignment * alignment, double * score,  double * res_fract_rank, int gap_treatment) {

    int *sorted_res;
    int pos, ctr, ctr2;
    int first, cvg_ctr;
    int almt_pos, gaps;
    int new_length = alignment->length;
    int usable_position, ref_s = 0, g;
    int *int_cvg;
    double *protein_score, prev_score; /* "score" refers to alignment positions */
    Group  *group, *ref_group;

    if (gap_treatment > 0) {
     	group = alignment->group+gap_treatment-1;
    }
    
    sorted_res    =    (int *) emalloc ( alignment->length*sizeof (int) );
    if (!sorted_res) return 1;

    if ( ! (int_cvg = emalloc(alignment->length*sizeof(int) ) ) ) return 1;
 
    if ( gap_treatment  <= -100 ) {
	g = -(gap_treatment+100);
	ref_group = alignment->group+ g;
	ref_s     = ref_group->group_member[0];
	gaps = 0;
	for (almt_pos=0; almt_pos < alignment->length; almt_pos ++ ) {
	    if ( score[almt_pos] == ALL_GAPS
		 || alignment->sequence[ref_s][almt_pos] == '.') gaps++;
	}
	
	new_length = alignment->length - gaps;
	
    } else if ( gap_treatment  == -2 ) {
	new_length = alignment->length;
	
    } else if ( gap_treatment == -1 ) {
	new_length = alignment->length - alignment->number_of_sunk_positions;
	
    } else if ( gap_treatment == 0 ) {
	gaps = 0;
	for (almt_pos=0; almt_pos < alignment->length; almt_pos ++ ) {
	    if ( score[almt_pos] == ALL_GAPS ) gaps++;
	}
	new_length = alignment->length - gaps;
	
    } else {
	ref_group = alignment->group+ gap_treatment-1;
	ref_s     = ref_group->group_member[0];
	gaps = 0;
	for (almt_pos=0; almt_pos < alignment->length; almt_pos ++ ) {
	    if (alignment->sequence[ref_s][almt_pos] == '.') gaps++;
	}
	new_length = alignment->length - gaps;
   }
 
    /*allocate */
    protein_score = (double *) emalloc ( new_length*sizeof (double) );
    if (!protein_score ) return 1;

     /* remove gapped positions from the score array */
    pos = 0;
    for (almt_pos=0; almt_pos < alignment->length; almt_pos ++ ) {


	usable_position = 1;
	if (gap_treatment <= -100) {
	    if ( score[almt_pos] == ALL_GAPS ||
		 alignment->sequence[ref_s][almt_pos] == '.')
		usable_position = 0;
	    
	} else if ( gap_treatment == -1  ) {
	    if ( alignment->sunk[almt_pos] )  usable_position = 0;
	    
	} else if (gap_treatment == 0) {
	    if ( score[almt_pos] == ALL_GAPS) usable_position = 0;
	    
	} else if (gap_treatment > 0) {
	    if (alignment->sequence[ref_s][almt_pos] == '.') usable_position = 0;
	    
	}
	
	if ( ! usable_position) continue;

	if ( pos >= new_length ) {
	    fprintf (stderr, "Error: pos ctr (%d) >= allocated length (%d) \n", pos, new_length);
	    exit (1);
	}
	protein_score[pos] = score[almt_pos];
	pos ++;
   }

    
    /* sort protein residues according to the new array */
    for (pos=0; pos < new_length; pos++) sorted_res[pos] = pos;
    array_qsort ( sorted_res, protein_score, new_length);

    /* turn the sorted array to coverage info */
    
    /* the lowest score in the game */
    prev_score = protein_score[ sorted_res[0] ];

    
    first   = 0;
    cvg_ctr = 0;
    int_cvg[cvg_ctr] = 0;
   
    for (ctr=0; ctr < alignment->length; ctr++) {
	res_fract_rank[ctr] = 1.0;
    }
    for (ctr=0; ctr < new_length; ctr++) {
	
	if ( protein_score[ sorted_res[ctr] ] <= prev_score ) {
	    int_cvg[cvg_ctr] ++;
	} else {
	    prev_score  = protein_score[ sorted_res[ctr] ];
	    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
		res_fract_rank[ sorted_res[ctr2] ] = (double) int_cvg[cvg_ctr]/new_length;
	    }
	    first = ctr;
	    cvg_ctr ++;
	    if ( cvg_ctr <  alignment->length) {
		int_cvg[cvg_ctr] =  int_cvg[cvg_ctr-1] + 1;
	    }
	}
    }
    for (ctr2=first; ctr2 <ctr; ctr2++ ) {
	res_fract_rank[ sorted_res[ctr2] ] =  (double) int_cvg[cvg_ctr]/new_length;
    }
    
 
    /* free */
    free (protein_score);
    free (sorted_res);
    free (int_cvg);
    
    return 0;
}
