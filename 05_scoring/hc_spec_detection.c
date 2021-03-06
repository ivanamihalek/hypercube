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


typedef enum {EUCLID, LINEAR} anon1;
typedef enum {NONE, ENTROPY, JS, ENTROPY_CORR,  MAX_CNSVTN_SCORE} anon2;
typedef enum {OVLP, FREQ_DIFF_COMPL, OVLP_CORR, MI_PAIR, JS_GROUP, MAX_OVERLAP_SCORE} anon3;


/***********************************************************************/
int evolve (double ** trans_prob, double *distr_in, double* distr_out);
int output_cube_coords (int no_pos, int no_groups, double **entropy, double *** overlap);
int scale_overlap (double *** overlap, int length, int no_groups, int **all_gaps,
		   double min_ovlp, double max_ovlp, int complement );

/* scoring functions */
int discriminants (int no_pos, int no_groups, int ** all_gaps,
		   double **variability, double *** overlap, int add_mode,  double *score);
int determinants (int no_pos, int no_groups, int ** all_gaps,
		  double **variability, double *** overlap, int add_mode,  int reference_group,  double *score);
int free_single_spec_distribution( double *** prob_matrix_series, int no_timesteps,
				   int no_of_aa_types, double ***single_spec_distr);
int model_evo_quantifiers(int no_pos, int no_groups, int ** all_gaps,
			  double ***freq, double ***freq2, 
			  int no_timesteps, int no_of_aa_types, double ***single_spec_distr,
			  double *avg_apex_time, double *max_apex_time, int **largest_free_evo_ovlp);

/***********************************************************************/
int  spec_detection ( Alignment * alignment, double *** prob_matrix_series, int no_timesteps,
		      double * stnry_freq, double **spec_score, char ** score_label, int *scores_assigned) {

    char label0, label1, label2;
    int pos, g, g1, g2, i, s, t, t2, aa, nongap_groups;
    int ancestral_aa, ancestral_aa2;
    int no_groups    = alignment->no_groups;
    int no_positions = alignment->length;
    int a_index;
    int score_ctr, cnsvtn_score, overlap_score, add_mode;
    int complement;
    int *nongap_col;
    int **largest_free_evo_ovlp;
    int **all_gaps;

    double norm, avg;
    double ovlp, diff;
    double correction, expected_entropy;
    double min_ovlp, max_ovlp;
    double min_entr  [no_groups];
    double max_entr  [no_groups];
    double min_entr_c[no_groups];
    double max_entr_c[no_groups];
 
    double *avg_apex_time;
    double *max_apex_time;
    double *expected_freq;
    double *JS_max;
    
    double **entropy;
    double **entropy_corr;
    double **JS_div;
    
    double ***overlap;
    double ***overlap_corr;
    double ***freq_diff_complement;
    double ***JS_group;
    double ***mi;
    double ***freq;
    double ***freq2;
    double ***single_spec_distr;

    Group *group;
 

    /* alloc */
    if ( ! (nongap_col =emalloc(no_groups*sizeof(int))) ) return 1;

    if ( ! (largest_free_evo_ovlp = intmatrix(no_positions,no_groups))) return 1;
    if ( ! (all_gaps  = intmatrix(no_positions,no_groups))) return 1;

    if ( ! (avg_apex_time = emalloc(no_groups*sizeof(double))) ) return 1;
    if ( ! (max_apex_time = emalloc(no_groups*sizeof(double))) ) return 1;
    if ( ! (JS_max        = emalloc(no_groups*sizeof(double))) ) return 1;
    if ( ! (expected_freq = emalloc(no_of_aa_types*sizeof(double))) ) return 1;

    if ( ! (entropy       = dmatrix(no_positions,no_groups))) return 1;
    if ( ! (entropy_corr  = dmatrix(no_positions,no_groups))) return 1;
    if ( ! (JS_div        = dmatrix(no_positions, no_groups))) return 1;
  
    if ( ! (overlap       = d3matrix(no_positions, no_groups, no_groups))) return 1;
    if ( ! (overlap_corr  = d3matrix(no_positions, no_groups, no_groups))) return 1;
    if ( ! (freq_diff_complement = d3matrix(no_positions, no_groups, no_groups))) return 1;
    if ( ! (JS_group      = d3matrix(no_positions, no_groups, no_groups))) return 1;
    if ( ! (mi            = d3matrix(no_positions, no_groups, no_groups))) return 1;
    if ( ! (freq          = d3matrix(no_positions, no_groups, no_of_aa_types))) return 1;
    if ( ! (freq2         = d3matrix(no_positions, no_groups, no_of_aa_types))) return 1;
    if ( ! (single_spec_distr = d3matrix(no_of_aa_types, no_timesteps, no_of_aa_types))) return 1;

    double **cons_ptr[MAX_CNSVTN_SCORE]     = {NULL, entropy, JS_div, entropy_corr};
    double***overlap_ptr[MAX_OVERLAP_SCORE] = {overlap, freq_diff_complement, overlap_corr, mi, JS_group};

    /****************************************************************/
    /****  counting the aa type frequencies                     *****/
    /****************************************************************/
    for (pos=0; pos< alignment->length; pos++) {
   
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
		if ( aa == 'X') continue;
		a_index = aa2index[aa];
       
		freq[pos][g][a_index] ++;
	    }
	    
	    memcpy    (freq2[pos][g], freq[pos][g], no_of_aa_types*sizeof(double));
	    normalize (freq [pos][g], no_of_aa_types, 1);
	    normalize (freq2[pos][g], no_of_aa_types, 2);

	    /* are we all gaps here ? */
	    a_index = aa2index['.'];
	    if ( freq[pos][g][a_index] > 0.3) {
		all_gaps[pos][g] = 1;
	    } else {
		all_gaps[pos][g] = 0;
	    }
	}
    }
    
    /****************************************************************/
    /****  'evolution' of a free system                         *****/
    /****************************************************************/
    if (options.exchangeability) {
	/* I'll need below:
	   single_spec_distr, ancestral aa and avg apex time            */

	/* the last aa type is '.' - we don't want to evolve that       */
	free_single_spec_distribution( prob_matrix_series,  no_timesteps,
				       no_of_aa_types-1, single_spec_distr);
	/* foreach group, for each position, find the majority type    */
	/* for the majority type find the time of the max overlap      */
	/* with the observed, and freely evolving distribution         */
	model_evo_quantifiers(alignment->length, no_groups, all_gaps, freq, freq2,
			      no_timesteps, no_of_aa_types-1, single_spec_distr,
			      avg_apex_time, max_apex_time, largest_free_evo_ovlp);
    }


    
    /****************************************************************/
    /****  entropy                                              *****/
    /****************************************************************/
    for (g=0; g<no_groups; g++) {
	min_entr[g] =  1000;
	max_entr[g] = -1000;
    }

    for (pos=0; pos< alignment->length; pos++) {
	for (g=0; g<no_groups; g++) {

	    if (all_gaps[pos][g] ) {
		entropy[pos][g] = 1.0;
		continue;
	    }
	    entropy[pos][g] = 0;
	    for (a_index=0; a_index<no_of_aa_types; a_index++) {
		if ( freq[pos][g][a_index] < 1.e-3) continue;
		entropy[pos][g] -= freq[pos][g][a_index]*log(freq[pos][g][a_index]);
	    }
	    //entropy[pos][g] /= log (20);
	    /* choose the scale for the entropy from the alignment itself: */
	    if (min_entr[g] > entropy[pos][g]) min_entr[g] = entropy[pos][g];
	    if (max_entr[g] < entropy[pos][g]) max_entr[g] = entropy[pos][g];
	}
    }
    /* entopy is rescaled below, after we have used it  to calculate the "corrected" entropy */

    
    /****************************************************************/
    /****  entropy corrected                                    *****/
    /****************************************************************/
    if (options.exchangeability) {
	for (g=0; g<no_groups; g++) {
	    min_entr_c[g] =  1000;
	    max_entr_c[g] = -1000;
	}
	for (pos=0; pos< alignment->length; pos++) {
	
	    for (g=0; g<no_groups; g++) {
		if ( all_gaps[pos][g]) continue;
		ancestral_aa = largest_free_evo_ovlp[pos][g];
		t = (int) avg_apex_time[g];
		memcpy (expected_freq, single_spec_distr[ancestral_aa][t], no_of_aa_types*sizeof(double));
	    
		/* we had it normalized to square above */
		normalize (expected_freq, no_of_aa_types, 1);
		expected_entropy = 0.0;
	    
		for (a_index=0; a_index<no_of_aa_types; a_index++) {
		
		    if ( expected_freq[a_index] >  1.e-3) 
			expected_entropy -= expected_freq[a_index]*log(expected_freq[a_index]);
		}

		if ( 1) {
		    entropy_corr[pos][g] = entropy[pos][g]-expected_entropy;
		} else {

		    if ( entropy[pos][g] && expected_entropy ) {
			correction = 1+log((10+entropy[pos][g])/(10+expected_entropy)) ;
		    } else {
			correction = 1;
		    }
		    entropy_corr[pos][g] = entropy[pos][g]*correction;

		}
		if (min_entr_c[g] > entropy_corr[pos][g]) min_entr_c[g] = entropy_corr[pos][g];
		if (max_entr_c[g] < entropy_corr[pos][g]) max_entr_c[g] = entropy_corr[pos][g];
	    }
	}
    
	for (pos=0; pos< alignment->length; pos++) {
	    for (g=0; g<no_groups; g++) {
		entropy[pos][g]      = (entropy[pos][g] - min_entr[g])/(max_entr[g] - min_entr[g]);
		entropy_corr[pos][g] = (entropy_corr[pos][g] - min_entr_c[g])/(max_entr_c[g] - min_entr_c[g]);
	    }
	}
    }
    
    /****************************************************************/
    /****  JS divergence                                        *****/
    /****************************************************************/
    if (options.compare ) {
	for (g=0; g<no_groups; g++) JS_max[g]      = -1000;

	for (pos=0; pos< alignment->length; pos++) {
	    for (g=0; g<no_groups; g++) {
		if ( all_gaps[pos][g]) continue;
       
		/*  JS: distance from the stationary distribution */
		JS_div[pos][g] = 0.0;
		for (a_index=0; a_index<no_of_aa_types; a_index++) {
		    if ( freq[pos][g][a_index]<1.e-6 || stnry_freq[a_index]<1.e-6) continue;
		    avg = (stnry_freq[a_index]+freq[pos][g][a_index])/2;
		    JS_div[pos][g] +=  freq[pos][g][a_index]*log(freq[pos][g][a_index]/avg);
		    JS_div[pos][g] +=  stnry_freq[a_index]*log(stnry_freq[a_index]/avg);
		}
		if (JS_max[g] < JS_div[pos][g]) JS_max[g] = JS_div[pos][g];
	    }
	
	    for (g=0; g<no_groups; g++) JS_div[pos][g] = 1-JS_div[pos][g]/JS_max[g];
	
	}
    }

    

    /****************************************************************/
    /****  direct overlaps between groups                       *****/
    /****************************************************************/
    min_ovlp =  1000;
    max_ovlp = -1000;
    for (pos=0; pos< alignment->length; pos++) {
	
	for (g=0; g<no_groups; g++) {
	    if (all_gaps[pos][g] ) continue;
       
	    for (g2=0; g2<no_groups; g2++) {
		if (all_gaps[pos][g2] ) continue;
		
		overlap[pos][g][g2] = 0.0;
		for (a_index=0; a_index<no_of_aa_types; a_index++) {
		    overlap[pos][g][g2] += freq2[pos][g][a_index]*freq2[pos][g2][a_index];
		}
		
		if (min_ovlp > overlap[pos][g][g2]) min_ovlp = overlap[pos][g][g2];
		if (max_ovlp < overlap[pos][g][g2]) max_ovlp = overlap[pos][g][g2];
	    }
	}
    }
    //scale_overlap (overlap, alignment->length, no_groups, all_gaps,  min_ovlp, max_ovlp, (complement=0) );
    
    /****************************************************************/
    /****  Jensen-Shannon divergence between groups             *****/
    /****************************************************************/
    if (options.compare ) {
	min_ovlp =  1000;
	max_ovlp = -1000;
	for (pos=0; pos< alignment->length; pos++) {
	
	    for (g=0; g<no_groups; g++) {
		if (all_gaps[pos][g] ) continue;
       
		for (g2=0; g2<no_groups; g2++) {
		    if (all_gaps[pos][g2] ) continue;
		
		    JS_group[pos][g][g2] = 0.0;
		    for (a_index=0; a_index<no_of_aa_types; a_index++) {
			avg = (freq[pos][g][a_index]+freq[pos][g2][a_index]);
			if ( avg < 1.0-3) continue;
			if ( freq[pos][g][a_index] ) {
			    JS_group[pos][g][g2] += freq[pos][g][a_index]*log(freq[pos][g][a_index] /avg);
			}
			if ( freq[pos][g2][a_index] ) {
			    JS_group[pos][g][g2] += freq[pos][g2][a_index]*log(freq[pos][g2][a_index] /avg);
			}
		    
		    }
		    if (min_ovlp > JS_group[pos][g][g2]) min_ovlp = JS_group[pos][g][g2];
		    if (max_ovlp < JS_group[pos][g][g2]) max_ovlp = JS_group[pos][g][g2];
		}
	    }
	}
	scale_overlap (JS_group, alignment->length, no_groups, all_gaps,  min_ovlp, max_ovlp, (complement=1) );
    }
    
    /****************************************************************/
    /****  overlap corrected                                    *****/
    /****************************************************************/
    if (options.exchangeability) {
	min_ovlp =  1000;
	max_ovlp = -1000;
    
	for (pos=0; pos< alignment->length; pos++) {
	
	    for (g=0; g<no_groups; g++) {
		if ( all_gaps[pos][g]) continue;
		ancestral_aa = largest_free_evo_ovlp[pos][g];
		t = (int) avg_apex_time[g];
		memcpy (expected_freq, single_spec_distr[ancestral_aa][t], no_of_aa_types*sizeof(double));
	    
		/* we had it normalized to square above */
		normalize (expected_freq, no_of_aa_types, 1);
		expected_entropy = 0.0;
	    
		for (g2=g+1; g2<no_groups; g2++) {
		    if ( all_gaps[pos][g2]) continue;
		    norm ++;

		    ancestral_aa2 = largest_free_evo_ovlp[pos][g2];
		    t2 = (int) avg_apex_time[g2];
       
		    ovlp = 0.0;
		    for (a_index=0; a_index<no_of_aa_types; a_index++) {
			ovlp += single_spec_distr[ancestral_aa][t][a_index]*
			    single_spec_distr[ancestral_aa2][t2][a_index];
		    }

		    if (0) {
			if (overlap[pos][g][g2] && ovlp) {
			    correction = 1+log(overlap[pos][g][g2]/ovlp) ;
			} else {
			    correction = 1.0;
			}
			overlap_corr[pos][g][g2] = overlap[pos][g][g2]*correction;
		    } else  {
			overlap_corr[pos][g][g2] = overlap[pos][g][g2]- ovlp;
		   
		    } 

		    
		    if (min_ovlp > overlap_corr[pos][g][g2]) min_ovlp = overlap_corr[pos][g][g2];
		    if (max_ovlp < overlap_corr[pos][g][g2]) max_ovlp = overlap_corr[pos][g][g2];

		}
	    }

	}

	scale_overlap (overlap_corr, alignment->length, no_groups, all_gaps,  min_ovlp, max_ovlp, (complement=0));
    }
     
    
    /****************************************************************/
    /****  freq difference between groups                       *****/
    /****************************************************************/
    if (options.compare ) {
	min_ovlp =  1000;
	max_ovlp = -1000;
	for (pos=0; pos< alignment->length; pos++) {
	    for (g=0; g<no_groups; g++) {
		if (all_gaps[pos][g] ) continue;
       
		for (g2=0; g2<no_groups; g2++) {
		    if (all_gaps[pos][g2] ) continue;
		
		    freq_diff_complement[pos][g][g2] = 0.0;
		    for (a_index=0; a_index<no_of_aa_types; a_index++) {
			diff = freq[pos][g][a_index]-freq[pos][g2][a_index];
			freq_diff_complement[pos][g][g2] += diff*diff; /*  Capra sticks in 0.5 here*/
		    }
		    if (min_ovlp > freq_diff_complement[pos][g][g2]) min_ovlp = freq_diff_complement[pos][g][g2];
		    if (max_ovlp < freq_diff_complement[pos][g][g2]) max_ovlp = freq_diff_complement[pos][g][g2];
		}
	    }
	}
	scale_overlap (freq_diff_complement, alignment->length, no_groups, all_gaps,  min_ovlp, max_ovlp, (complement=1));

    }
    
    /****************************************************************/
    /****  pairwise mutual info                                *****/
    /****************************************************************/
    if (options.compare ) {
	min_ovlp =  1000;
	max_ovlp = -1000;
	for (pos=0; pos< alignment->length; pos++) {

	    int element;
	    int pair[2];
	    double mutual_info, fract_group_size;
	    double all_freq[no_of_aa_types];
	    double freq_of_aa_group_pair[no_groups][no_of_aa_types];
	

	    for (g1=0; g1<no_groups; g1++) {
		if (all_gaps[pos][g1] ) continue;
		pair[0] = g1;
	    
		for (a_index=0; a_index<no_of_aa_types; a_index++) {
		    freq_of_aa_group_pair[g1][a_index] = 0.0;
		}
   
		for (g2=g1+1; g2<no_groups; g2++) {
		    if (all_gaps[pos][g2] ) continue;
		    pair[1] = g2;

		    /* overall freq distribution in the two groups */
		    for (a_index=0; a_index<no_of_aa_types; a_index++) all_freq[a_index] = 0;
		    for (element=0; element<2; element++) {
			g     = pair[element];
			group = alignment->group+g;
			for (i=0; i< group->no_members; i++) {
			    s  = group->group_member[i];
			    aa = (int) alignment->sequence[s][pos];
			    /* skip x when counting frequencies  */
			    if ( aa == 'X') continue;
			    a_index = aa2index[aa];
			    all_freq[a_index] ++;
			}
		    }
		    normalize (all_freq, no_of_aa_types, 1);

		    /* distribution of aa-group pairing */
		    int  number_of_seqs_in_the_two_groups = 0;
		    for (element=0; element<2; element++) {
			g     = pair[element];
			group = alignment->group+g;
			number_of_seqs_in_the_two_groups += group->no_members;
		    }
		    for (element=0; element<2; element++) {
			g = pair[element];
			for (a_index=0; a_index<no_of_aa_types; a_index++) {
			    freq_of_aa_group_pair[g][a_index] = 0.0;
			}
			group = alignment->group+g;
			for (i=0; i< group->no_members; i++) {
			    s  = group->group_member[i];
			    aa = (int) alignment->sequence[s][pos];
			    /* skip x when counting frequencies  */
			    if ( aa == 'X') continue;
			    a_index = aa2index[aa];
			    freq_of_aa_group_pair[g][a_index] ++;
			}
			for (a_index=0; a_index<no_of_aa_types; a_index++) {
			    freq_of_aa_group_pair[g][a_index] /=  number_of_seqs_in_the_two_groups;
			}
		    }

	
		    /* mi calculation */
		    mutual_info = 0.0;
		
		    for (element=0; element<2; element++) {
			g     = pair[element];
			group = alignment->group+g;
			fract_group_size = (double)group->no_members/number_of_seqs_in_the_two_groups;
	    
			for (a_index=0; a_index<no_of_aa_types; a_index++) {
			    if (!freq_of_aa_group_pair[g][a_index] || !all_freq[a_index] )  continue;
			    mutual_info += freq_of_aa_group_pair[g][a_index]
				*log(freq_of_aa_group_pair[g][a_index]/(all_freq[a_index]*fract_group_size));
			}
		    }
		    mi[pos][g1][g2] =  mi[pos][g2][g1] = mutual_info;
		    if (min_ovlp > mutual_info) min_ovlp = mutual_info;
		    if (max_ovlp < mutual_info) max_ovlp = mutual_info;

		}
	    }  
	}  
	scale_overlap (mi, alignment->length, no_groups, all_gaps,  min_ovlp, max_ovlp, (complement=1) );
    }

    
    /****************************************************************/
    /****  output cube coords                                   *****/
    /****************************************************************/
    if ( 0 ) output_cube_coords (alignment->length, no_groups, entropy, overlap);


    /****************************************************************/
    /****************************************************************/
    /****************************************************************/
    /****************************************************************/
    /****************************************************************/
    /****************************************************************/
    score_ctr = 0;

    if (options.compare ) {
	/****************************************************************/
	/****  calculate various scoring combos                    *****/
	/****************************************************************/
    
	for (cnsvtn_score=0; cnsvtn_score<MAX_CNSVTN_SCORE; cnsvtn_score++) {
 	
	    if (cnsvtn_score==NONE ) {label0='0';}
	    else if (cnsvtn_score==ENTROPY) {label0='e';}
	    else if (cnsvtn_score==JS) {label0='j';}
	    else if (cnsvtn_score==ENTROPY_CORR) {label0='r';}
	    
	    for (overlap_score=0; overlap_score<MAX_OVERLAP_SCORE; overlap_score++) {

		if (overlap_score==OVLP) {label1='o';}
		else if (overlap_score==FREQ_DIFF_COMPL) {label1='f';}
		else if (overlap_score==OVLP_CORR) {label1='r';}
		else if (overlap_score==MI_PAIR)   {label1='m';}
		else if (overlap_score==JS_GROUP)  {label1='j';}
	    
		score_label[score_ctr][1] = label1;
	    
		for (add_mode=EUCLID; add_mode<=LINEAR; add_mode++) {

		    if (add_mode==EUCLID) {label2='e';}
		    else if (add_mode==LINEAR) {label2='l';}
		    score_label[score_ctr][2] = label2;

		
		    if (score_ctr >= MAX_NO_SCORES) {
			fprintf (stderr, "underdimensioned score array (increase MAX_NO_SCORES)\n");
			exit (1);
		    }
		
		    sprintf (score_label[score_ctr], "%c%c%c%s", label0, label1, label2, "_ds");
		    discriminants (alignment->length, no_groups, all_gaps,
				   cons_ptr[cnsvtn_score], overlap_ptr[overlap_score], add_mode, spec_score[score_ctr]);
		    score_ctr ++;
		
		    if (score_ctr >= MAX_NO_SCORES) {
			fprintf (stderr, "underdimensioned score array (increase MAX_NO_SCORES)\n");
			exit (1);
		    }
		    sprintf (score_label[score_ctr], "%c%c%c%s", label0, label1, label2, "_dt");

		    determinants  (alignment->length, no_groups, all_gaps,
				   cons_ptr[cnsvtn_score], overlap_ptr[overlap_score],  add_mode, 0, spec_score[score_ctr]);
		    score_ctr ++;
		
		}
	    }
	}
    
	/****************************************************************/
	/****  mutual information                                   *****/
	/****************************************************************/
	if (score_ctr >= MAX_NO_SCORES) {
	    fprintf (stderr, "underdimensioned score array (increase MAX_NO_SCORES)\n");
	    exit (1);
	}
   
	sprintf (score_label[score_ctr], " MI ");
    
	for (pos=0; pos< alignment->length; pos++) {
	
	    double mutual_info, fract_group_size;
	    double all_freq[no_of_aa_types];
	    double freq_of_aa_group_pair[no_groups][no_of_aa_types];
	

	    for (a_index=0; a_index<no_of_aa_types; a_index++) all_freq[a_index] = 0;
	    for (s=0; s<alignment->number_of_seqs; s++) {
		aa = (int) alignment->sequence[s][pos];
		/* skip x when counting frequencies  */
		if ( aa == 'X') continue;
		a_index = aa2index[aa];
       
		all_freq[a_index] ++;
	    }
	    normalize (all_freq, no_of_aa_types, 1);
	
	    for (g=0; g<no_groups; g++) {
		for (a_index=0; a_index<no_of_aa_types; a_index++) {
		    freq_of_aa_group_pair[g][a_index] = 0.0;
		}
	    }

	    for (g=0; g<no_groups; g++) {
		group = alignment->group+g;
		for (i=0; i< group->no_members; i++) {
		    s  = group->group_member[i];
		    aa = (int) alignment->sequence[s][pos];
		    /* skip x when counting frequencies  */
		    if ( aa == 'X') continue;
		    a_index = aa2index[aa];
		    freq_of_aa_group_pair[g][a_index] ++;
		}
	    }

	    for (g=0; g<no_groups; g++) {
		for (a_index=0; a_index<no_of_aa_types; a_index++) {
		    freq_of_aa_group_pair[g][a_index] /= alignment->number_of_seqs;
		}
	    }

	
	    mutual_info = 0.0;
	    for (g=0; g<no_groups; g++) {

		group = alignment->group+g;
		fract_group_size = (double)group->no_members/alignment->number_of_seqs;
	    
		for (a_index=0; a_index<no_of_aa_types; a_index++) {
		    if (!freq_of_aa_group_pair[g][a_index] || !all_freq[a_index] )  continue;
		    /* neg sign to conform to the same order as the other methods */ 
		    mutual_info -= freq_of_aa_group_pair[g][a_index]
			*log(freq_of_aa_group_pair[g][a_index]/(all_freq[a_index]*fract_group_size));
		}
    
	    }
	
	    spec_score[score_ctr][pos] = mutual_info;
	}
	score_ctr ++;
	
    } else {
	/**********************************************************************/
	/****  we are not doing the across-the board comparison,          *****/
	/****  we just need a working version of dets and discriminants   *****/
	/**********************************************************************/
	
	add_mode = LINEAR;
	
	if (options.exchangeability) {
	    cnsvtn_score  = ENTROPY_CORR;
	    overlap_score = OVLP_CORR;
	} else {
	    cnsvtn_score  = ENTROPY;
	    overlap_score = OVLP;
	}
	discriminants (alignment->length, no_groups, all_gaps,
		       cons_ptr[cnsvtn_score], overlap_ptr[overlap_score],
		       add_mode, spec_score[score_ctr]);
	score_ctr ++;
		
	for (g=0; g<no_groups; g++) {
	    determinants  (alignment->length, no_groups, all_gaps,
			   cons_ptr[cnsvtn_score], overlap_ptr[overlap_score],
			   add_mode, g, spec_score[score_ctr]);
	    score_ctr ++;
	}

    }
    
   
    *scores_assigned = score_ctr;

    
    /****************************************************************/
    /****  for all scores, mark the gapped positions            *****/
    /****************************************************************/
    for (pos=0; pos< alignment->length; pos++) {

	nongap_groups = 0;
	for (g=0; g<no_groups; g++) {
	    if ( all_gaps[pos][g]) continue;
	    nongap_groups ++;
	}

	if (  nongap_groups <=1) {
	    for (score_ctr = 0; score_ctr< *scores_assigned; score_ctr++) {
		spec_score[score_ctr][pos] = ALL_GAPS;
	    }
	    continue;
	}
    }
    score_ctr = 0;


    
    /****************************************************************/
    /*     cleanup                                                  */
    /****************************************************************/
    free (nongap_col);
    
    free_imatrix (largest_free_evo_ovlp);
    free_imatrix (all_gaps);

    free(avg_apex_time);
    free(max_apex_time);
    free(expected_freq);
    free(JS_max);

    free_dmatrix (entropy);
    free_dmatrix (entropy_corr);
    free_dmatrix (JS_div);

    free_d3matrix(overlap);
    free_d3matrix(overlap_corr);
    free_d3matrix(freq_diff_complement);
    free_d3matrix(JS_group);
    free_d3matrix(mi);
    free_d3matrix(freq);
    free_d3matrix(freq2);
    free_d3matrix(single_spec_distr);

    
    return 0;
}


/******************************************************/
/******************************************************/
/******************************************************/
int discriminants (int no_pos, int no_groups, int ** all_gaps,
		   double **variability, double *** overlap, int add_mode,  double *score) {
    
    int pos, g, g2;
    int eff_no_pairs, eff_no_groups;
    double avg_ovlp, avg_variability;
    double discr;
    
    for (pos=0; pos< no_pos; pos++) {

	discr = 0.0;
	eff_no_pairs = eff_no_groups   = 0;
	avg_ovlp     = avg_variability = 0.0;
	
	for (g=0; g<no_groups; g++) {
	    if ( all_gaps[pos][g]) continue;
       
	    for (g2=g+1; g2<no_groups; g2++) {
		if ( all_gaps[pos][g2]) continue;
		if (add_mode == EUCLID) {
		    avg_ovlp +=  overlap[pos][g][g2]*overlap[pos][g][g2];
		} else {
		    avg_ovlp +=  overlap[pos][g][g2];
		}
		eff_no_pairs ++;
	    }
	    if ( variability ) {
		if ( add_mode == EUCLID) {
		    avg_variability += variability[pos][g]*variability[pos][g];
		} else {
		    avg_variability += variability[pos][g];
		}
	    }
	    eff_no_groups ++;
	}
	avg_ovlp /= eff_no_pairs;
	avg_variability /= eff_no_groups;
	
	discr = 2*avg_ovlp + avg_variability;
	if ( add_mode == EUCLID) discr = sqrt(discr);
	score[pos] =  discr;
    }

    return 0;
}

/******************************************************/
int determinants (int no_pos, int no_groups, int ** all_gaps,
		  double **variability, double *** overlap, int add_mode,  int reference_group,  double *score) {
    
    int pos, g2;
    int eff_no_pairs, eff_no_groups;
    double avg_ovlp, avg_variability;
    double det;
    
    for (pos=0; pos< no_pos; pos++) {

	if ( all_gaps[pos][reference_group]) continue;
	det = 0.0;
	eff_no_pairs = eff_no_groups = 0;
	avg_ovlp = avg_variability = 0.0;
	
	for (g2=0; g2<no_groups; g2++) {
	    if (g2== reference_group) continue;
	    if (all_gaps[pos][g2]) continue;
	    
	    if (add_mode == EUCLID) {
		avg_ovlp +=  overlap[pos][reference_group][g2]*overlap[pos][reference_group][g2];
	    } else {
		avg_ovlp +=  overlap[pos][reference_group][g2];
	    }
	    eff_no_pairs ++;

	}
	if ( variability ) {
	    if (add_mode == EUCLID) {
		avg_variability += variability[pos][reference_group]*variability[pos][reference_group];
	    } else {
		avg_variability += variability[pos][reference_group];
	    }
	}
	eff_no_groups ++;
	
	avg_ovlp /= eff_no_pairs;
	avg_variability /= eff_no_groups;
	
	det = avg_ovlp + avg_variability;
	if ( add_mode == EUCLID) det = sqrt(det);
	score[pos] =  det;


    }
   
    return 0;
}




/******************************************************/
/******************************************************/

int evolve (double ** trans_prob, double *distr_in, double* distr_out) {
    int a_index, b_index;
    for (a_index=0; a_index<no_of_aa_types; a_index++) {
	distr_out[a_index] = 0.0;
	for (b_index=0; b_index<no_of_aa_types-1; b_index++) {
	    distr_out[a_index] += trans_prob[a_index][b_index]*distr_in[b_index];
	}
    }
    return 0;
   
}


/******************************************************/
/******************************************************/
int normalize (double *array, int size, int power) {

    int i, p;

    double norm = 0, aux;

    for (i=0; i<size; i++) {
	aux = 1;
	for (p=0; p<power; p++) {
	    aux *= array[i];
	}
	norm += aux;
    }

    if (power > 2) {
	fprintf (stderr, "%d-norm not implemented\n", power);
	exit (1);
    } else if ( power== 2 ) {
	norm = sqrt (norm);
    }

    if ( norm < 1.e-6) return 0;
    
    for (i=0; i<size; i++) {
	array[i] /= norm;
    }
   
    return 0;
}

/******************************************************/
/******************************************************/
int scale_overlap (double *** overlap, int length, int no_groups, int **all_gaps,
		   double min_ovlp, double max_ovlp, int complement ) {
    int pos, g, g2;
    
    for (pos=0; pos< length; pos++) {
	for (g=0; g<no_groups; g++) {
	    if ( all_gaps[pos][g]) continue;
	    for (g2=g+1; g2<no_groups; g2++) {
		if ( all_gaps[pos][g2]) continue;
		overlap[pos][g][g2] = (overlap[pos][g][g2]-min_ovlp)/(max_ovlp-min_ovlp);
		if ( complement )  overlap[pos][g][g2] = 1.0 - overlap[pos][g][g2];
		overlap[pos][g2][g] = overlap[pos][g][g2];
	    }
	}
    }
    return 0;
}

/******************************************************/
/******************************************************/
int free_single_spec_distribution( double *** prob_matrix_series, int no_timesteps,
				   int no_of_aa_types, double ***single_spec_distr){

    int a_index, b_index, t;
    
    for (a_index=0; a_index<no_of_aa_types; a_index++) {

	if ( amino_acid_order[a_index] == '.') continue;
	
	
	for (b_index=0; b_index<no_of_aa_types; b_index++) {
	    single_spec_distr[a_index][0][b_index] = 0.1;/* bonferoni ") */
	    //single_spec_distr[a_index][0][b_index] = 0.0;
	}
	single_spec_distr[a_index][0][a_index] = 1.0;
	normalize (single_spec_distr[a_index][0], no_of_aa_types, 2);

	for (t=1; t<no_timesteps; t+=1) {
	    memset (single_spec_distr[a_index][t],  0, no_of_aa_types*sizeof(double));
	    evolve (prob_matrix_series[t], single_spec_distr[a_index][0], single_spec_distr[a_index][t]);
	    /* normalize to square, bcs that's what we'll be needing before */
	    normalize (single_spec_distr[a_index][t], no_of_aa_types, 2);
	}
    }

    
    return 0;
}

/******************************************************/
/******************************************************/
int output_cube_coords (int no_pos, int no_groups, double **entropy, double *** overlap) {
    
    int pos, g, g2;
    FILE * fptr;
    
    fptr = efopen ("tmp.no_sim.cube","w");
    if ( ! fptr) {
	fprintf (stderr, "error opening tmp.no_sim.cube.\n");
	exit (1);
    }
    for (pos=0; pos< no_pos; pos++) {
	fprintf (fptr, "%4d  ", pos+1);
	for (g=0; g<no_groups; g++) {
	    fprintf (fptr, " %8.3lf ", entropy[pos][g]);
	}
	for (g=0; g<no_groups; g++) {
	    for (g2=g+1; g2<no_groups; g2++) {
		fprintf (fptr, " %8.3lf ", overlap[pos][g][g2]);
	    }
	}
	fprintf (fptr, "\n");
    }
    fclose (fptr);

    return 0;
}

/******************************************************/
/******************************************************/
int model_evo_quantifiers(int no_pos, int no_groups, int ** all_gaps,
			  double ***freq, double ***freq2, 
			  int no_timesteps, int no_of_aa_types, double ***single_spec_distr,
			  double *avg_apex_time, double *max_apex_time, int **largest_free_evo_ovlp){

    int g, pos, a_index, j, t;
    int ancestral_aa;
    int nongap_col[no_groups];
    int *maj_a_index;
    int **apex_time;
    double ovlp,** apex_ovlp;
    double *maj_fraction;
    
    if ( ! (maj_a_index=emalloc(no_of_aa_types*sizeof(int))) ) return 1;
    if ( ! (apex_time = intmatrix(no_pos,no_groups))) return 1;
    if ( ! (apex_ovlp = dmatrix(no_pos,no_groups))) return 1;
    if ( ! (maj_fraction  = emalloc(no_of_aa_types*sizeof(double))) ) return 1;
    /* foreach group, for each position, find the majority type    */
    /* for the majority type find the time of the max overlap      */
    /* with the observed, and freely evolving distribution         */
    for (g=0; g<no_groups; g++) {
	
	avg_apex_time[g] = 0;
	max_apex_time[g] = 0;
	nongap_col[g] = 0;
	for (pos=0; pos< no_pos; pos++) {
	    if (all_gaps[pos][g] ) continue;

	    memset (maj_fraction, 0, no_of_aa_types*sizeof(double) );
	    memset ( maj_a_index, 0, no_of_aa_types*sizeof(int) );

	    /* sort the types by frequency */
	    for (a_index=0; a_index<no_of_aa_types; a_index++) {
		if ( amino_acid_order[a_index] == '.' ) continue;
		if (!freq[pos][g][a_index]) continue;
		j = 0;
		while (j<no_of_aa_types && maj_fraction[j]>freq[pos][g][a_index]) j++;
		if ( j<no_of_aa_types ) {
		    /* move*/
		    memmove (maj_fraction+j+1, maj_fraction+j, (no_of_aa_types-1-j)*sizeof(double));
		    maj_fraction[j] = freq[pos][g][a_index];
           
		    memmove (maj_a_index+j+1, maj_a_index+j, (no_of_aa_types-1-j)*sizeof(int));
		    maj_a_index[j]       = a_index;
		}
	    }
	    
	    /* at which time does the overlap have max ? */
	    apex_time[pos][g] = 0;
	    apex_ovlp[pos][g] = -1.0;
	    j = 0;
	    while (j<no_of_aa_types && maj_fraction[j]>= 0.9*maj_fraction[0] ) {

		double apex_freq[no_of_aa_types];
       
		/* the distribution consisits of maj_a_index[j] only */
		ancestral_aa = maj_a_index[j];
       
		for (t=0; t<no_timesteps; t+=1) {
		    /* calculate the overlap btw freely evolving my freq vector */
		    ovlp = 0;
		    for (a_index=0; a_index<no_of_aa_types; a_index++) {
			ovlp += single_spec_distr[ancestral_aa][t][a_index]*freq2[pos][g][a_index];
		    }
		    if ( apex_ovlp[pos][g] < ovlp ) {
			apex_ovlp[pos][g] = ovlp;
			apex_time[pos][g] = t;
			memcpy (apex_freq, single_spec_distr[ancestral_aa][t], no_of_aa_types*sizeof(double));
			largest_free_evo_ovlp[pos][g] = ancestral_aa;
		    }
		} /* on to the next time point */

		/* on to the next maj_fraction aa, if there is a tie */
		j++;
       
	    }

	    avg_apex_time[g] += apex_time[pos][g];
	    if (max_apex_time[g] < apex_time[pos][g] )  max_apex_time[g] = apex_time[pos][g];
	    nongap_col[g] ++;
	    
	}
    }
    
    for (g=0; g<no_groups; g++) {
	avg_apex_time[g] /= nongap_col[g];
    }
     
    free (maj_fraction);
    free (maj_a_index);
    free_imatrix (apex_time);
    free_dmatrix (apex_ovlp);
 
    return 0;
    
}




# if 0
    /* junkyard */
    ancestral_aa = 10; /* that's leucine */
    for (t=0; t<no_timesteps; t+=1) {
	for (a_index=0; a_index<no_of_aa_types; a_index++) {
	    memcpy (expected_freq, single_spec_distr[ancestral_aa][t], no_of_aa_types*sizeof(double));
	    
	    /* we had it normalized to square above */
	    normalize (expected_freq, no_of_aa_types, 1);
	    printf (" %4d   %3d   %8.4f\n", t, a_index, expected_freq[a_index]);
	}
    }
	    
    exit (1);




 	g = 0; g2 = 1;
	for (pos=0; pos< alignment->length; pos++) {
	    if ( all_gaps[pos][g2] ||  all_gaps[pos][g]) continue;
	    printf ( "%4d  %d  %d     %8.3lf   %8.3lf       %8.3lf   %12.8lf      %8.3lf   %8.3lf  \n",
		     pos+1, g, g2, entropy[pos][g], entropy[pos][g2], 
		     overlap[pos][g][g2], overlap_corr[pos][g][g2], min_ovlp, max_ovlp);
	}	    



	if ( reference_group==1 && pos==84 ) {
	    printf (" %8.4lf   %8.4lf   %8.4lf   %8.4lf  \n",
		    avg_variability, overlap[pos][reference_group][0], overlap[pos][0][reference_group], det );
	}

# endif
