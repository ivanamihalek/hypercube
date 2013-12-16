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
#include "hypercube.h"


int  comparative_score_output (Alignment *alignment, double ** score, char**score_label, int no_scores){

    int almt_pos, score_ctr;
    int g, pos, s;
 
    char filename[BUFFLEN];
    char aa = '\0';
    double cvg;
    double ** res_fract_rank = NULL;
    FILE * fptr;
    Group* group;
    
    sprintf (filename, "%s.comp_score", options.outname);

    fptr = efopen (filename, "w");
    if ( !fptr) return 1;

    /* the method names */
    /* find residue ranks, unless we are explicitly 
       instructed to output raw scores */
    if ( !options.raw_score) {
	
	if ( !(res_fract_rank=dmatrix(no_scores, alignment->length) ) ) return 1;
	
	for ( score_ctr=0; score_ctr<no_scores; score_ctr++) {
	    for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {
		if ( score[0][almt_pos] == ALL_GAPS )  {
		    score[score_ctr][almt_pos] = ALL_GAPS;
		} 
	    }
	    coverage (alignment, score[score_ctr], res_fract_rank[score_ctr], 0);
	}

    }

    
    
    fprintf (fptr, "%% %4s", "almt");
    for (g=0; g<alignment->no_groups; g++) {
	group = alignment->group+g;
	s  = group->group_member[0];
	fprintf ( fptr, " %10s", alignment->name[s]);
    }
    fprintf (fptr, "%8s", "gaps");
    
    for ( score_ctr=0; score_ctr<no_scores; score_ctr++) {
	fprintf (fptr, " %8s ", score_label[score_ctr]);
    }
    fprintf (fptr, "\n");
    
    pos = -1;
    for (almt_pos = 0; almt_pos < alignment->length; almt_pos++) {

  	fprintf (fptr, "%6d", almt_pos+1);
	/*  representative seqs */
	for (g=0; g<alignment->no_groups; g++) {
	    group = alignment->group+g;
	    s  = group->group_member[0];
	    aa = alignment->sequence[s][almt_pos];
	    fprintf (fptr, "%3c ", aa);
	}
	fprintf (fptr, "%5.2lf", (double)alignment->column_gaps[almt_pos]/alignment->number_of_seqs);
	
	if ( score[0][almt_pos] != ALL_GAPS )  pos++;
	/* scores */
	for ( score_ctr=0; score_ctr<no_scores; score_ctr++) {
	    if ( score[0][almt_pos] == ALL_GAPS ) { /* should be the same for all scores */
		//fprintf (fptr, "%10s%10s", "-", "-");
		//fprintf (fptr, "%10s", "-");
		fprintf (fptr, " %4.2lf ", 2.0);
	    } else {
		//fprintf (fptr, "%8.2lf%8.2lf", score[score_ctr][almt_pos], reliability[score_ctr][almt_pos]);
		
		if (options.raw_score) {
		     fprintf (fptr, " %8.2lf ", score[score_ctr][almt_pos]);
		} else {
		    if (alignment->sunk[almt_pos]) {
			cvg = 1;
		    } else {
			cvg = res_fract_rank[score_ctr][pos];
		    }
		    fprintf (fptr, " %4.2lf ", cvg);
		}
	    }
	}
	fprintf (fptr, "\n");
    }

    fclose (fptr);

    if ( !options.raw_score) free_dmatrix (res_fract_rank);
    
    return 0;

}

