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


/******************************************************************************/
int conservation_scoring (Options *options, Alignment * alignment, 
			  int * similar_to, double ** score, double *** in_group_score) {

    int retval;
    Tree tree;
    
    int build_tree (Options * options, Alignment * alignment, Tree * tree);
    int entropy ( Alignment * alignment, int * similar_to, double *score);
    int hybrid (Alignment * alignment, Tree * tree,
		 int * similar_to, int normalize,  double *score);
    int in_group_entropy ( Alignment * alignment, int * similar_to, double **score);
			   
    /* build the seq similarity tree */
    //printf ("building the tree: ...");  fflush (stdout);
    memset (&tree, 0, sizeof(Tree));  
    retval  = build_tree(options, alignment, &tree);
    if (retval) return retval;

    hybrid  (alignment, &tree, NULL, 0, score[0]);
    entropy (alignment, NULL, score[1]);
 
    /* per-group scoring */
    in_group_entropy (alignment, NULL, in_group_score[0]);
    /* if I ever decide to add other cons methods, sub as rvet */
    /* it should be something like */
    /* in_group_rvET (alignment, NULL, in_group_score[1]); */
   
    
    /* sink the gapped positions to the bottom, if requested */
    if ( options->max_gaps ) {
	int method_ctr;
	int ctr, group;
	double max_score = -1;

	alignment->number_of_sunk_positions = 0;
	for ( ctr=0; ctr < alignment->length; ctr++ ) {
	    if ((double)alignment->column_gaps[ctr]/alignment->number_of_seqs > options->max_gaps ) {
		alignment->sunk[ctr] = 1;
		alignment->number_of_sunk_positions++;
	    }
	}
	for (method_ctr=0; method_ctr<=1; method_ctr++) {

	    max_score = -1;
	    
	    for ( ctr=0; ctr < alignment->length; ctr++ ) {
		if ( max_score < score[method_ctr][ctr] ) max_score = score[method_ctr][ctr];
	    }
	    for ( ctr=0; ctr < alignment->length; ctr++ ) {
		if (alignment->sunk[ctr] ) {
		    score[method_ctr][ctr] = max_score;
		}
	    }
	}
	/* in group entropy */
	for (group=0;  group< alignment->no_groups; group++) {

	    max_score = -1;
	    
	    for ( ctr=0; ctr < alignment->length; ctr++ ) {
		if (max_score < in_group_score[0][group][ctr] )
		    max_score = in_group_score[0][group][ctr];
	    }
	    
	    for ( ctr=0; ctr < alignment->length; ctr++ ) {
		if ((double)alignment->group[group].gaps[ctr]/alignment->group[group].no_members
		    > options->max_gaps ) {
		    in_group_score[0][group][ctr] = max_score;
		}
	    }
	    
	}
    }

    free_node_matrix (tree.group_root);
    free (tree.leaf);
    return 0;
}

/******************************************************************************/
int entropy ( Alignment * alignment, int * similar_to, double *score){

    int col, seq, ctr, aa;
    int freq[ASCII], norm; /* ASCII == 128,  ASCII size  */
    double p, entropy;
    
    for (col = 0; col < alignment->length; col++){
	/* find frequencies */
	memset (freq, 0, ASCII*sizeof(int));
	norm = 0;
	for (seq=0; seq< alignment->number_of_seqs; seq++){
	    aa = (int) alignment->sequence[seq][col];
	    if ( aa == 'X' || aa == 'x' ) continue;
	    if (similar_to) aa=similar_to[aa];
	    freq [aa]++;
	    norm ++;
	}
	/* find entropy */
	entropy = 0.0;
	for ( ctr=0; ctr < ASCII; ctr++) {
	    if ( freq[ctr] ) {
		p = (double)freq[ctr]/norm;
		entropy -= p*log(p);
	    }
	}
	score[col] = entropy;
    }

    return 0;
}

/******************************************************************************/
int in_group_entropy ( Alignment * alignment, int * similar_to, double **score){

    int group, col, seq, ctr, aa;
    int **freq, *norm; /* ASCII == 128,  ASCII size  */
    double p, entropy;

    if ( ! (freq=intmatrix(alignment->no_groups, ASCII)) ) return 1;
    if ( ! (norm=emalloc(alignment->no_groups*sizeof(int))) ) return 1;
    
    for (col = 0; col < alignment->length; col++){

	/* find frequencies */
	for (group=0;  group< alignment->no_groups; group++) {
	    memset (freq[group], 0, ASCII*sizeof(int));
	    norm[group] = 0;
	}
	
	for (seq=0; seq< alignment->number_of_seqs; seq++){
	    aa = (int) alignment->sequence[seq][col];
	    if ( aa == 'X' || aa == 'x' ) continue;
	    if (similar_to) aa=similar_to[aa];
	    group = alignment->belongs_to_group[seq];
 	    freq [group][aa]++;
	    norm [group] ++;
	}
	/* find entropy */
	for (group=0;  group< alignment->no_groups; group++) {
	    entropy = 0.0;
	    for ( ctr=0; ctr < ASCII; ctr++) {
		if ( freq[group][ctr] ) {
		    p = (double)freq[group][ctr]/norm[group];
		    entropy -= p*log(p);
		}
	    }
	    score[group][col] = entropy;
	}
    }

    free_imatrix (freq);
    free (norm);
    
    return 0;
}

/******************************************************************************/
int hybrid (Alignment * alignment, Tree * tree, int * similar_to,
	    int normalize,  double *score){
    
    int pos, no_seqs = alignment->number_of_seqs;
    int rank, g;
    double subsum, rho;
    int entropy_recursive (Node *node,  int pos, int * similar_to, int normalize);
    
    /* calculate the score */
    for ( pos=0; pos < alignment->length; pos++) {
	entropy_recursive(tree->root, pos, similar_to, normalize);
	rho = 1.0;
	for ( rank = 1; rank < no_seqs; rank++ ) {
	    subsum = 0.0;
	    for ( g=0; g<rank; g++) {
		subsum  +=  tree->group_root[rank][g]->entropy;
	    }
	    rho += subsum/rank;
	}
 	score[pos] =  rho;
    }

    return 0;
    
}

/******************************************************************************/
int entropy_recursive (Node *node,  int pos, int * similar_to, int normalize) {

    if ( node->type == LEAF ) {
	int  aa;
	memset (node->population, 0, ASCII*sizeof(int));
	aa = (int) node->seq[pos];
	if ( similar_to) aa = similar_to[aa];
	node->population[ aa ] = 1;
	node->entropy = 0.0;
    } else {
	int ctr, norm;
	int no_of_leaves = node->number_of_leaves;
	int no_of_types = 0;
	double fr;
	entropy_recursive (node->left, pos, similar_to, normalize);
	entropy_recursive (node->right, pos, similar_to, normalize);
	/* for each aa type: population = population left + pop right */
	node->entropy = 0.0;
	if ( normalize ) {
	    ctr = '.';
	    node->population[ctr] = node->left->population[ctr]
		+ node->right->population[ctr];
	    no_of_leaves -= node->population[ctr];
	}
	norm = 0;
	for (ctr=0; ctr < ASCII; ctr++ ) {
	    if (normalize && ctr == '.') continue;
	    if ( ctr=='X' || ctr == 'x') continue;
	    node->population[ctr] = node->left->population[ctr]
		+ node->right->population[ctr];
	    norm += node->population[ctr];
	}
	for (ctr=0; ctr < ASCII; ctr++ ) {
	    /* frequency is pop/number of seqs */
	    if ( node->population[ctr] ) {
		if ( normalize ) no_of_types++;
		fr =  (double)node->population[ctr]/norm;
		node->entropy -= fr*log(fr);
	    }
	}
	if (normalize && no_of_types>1)
	    node->entropy /= log( no_of_types) ;
   }

    return 0;
}

/******************************************************************************/
# if 0
if ( col==80 || col==81) {
    printf ("\n col %d group %d \n", col, group);
    printf ("entr  %8.2lf \n", entropy);
}
# endif
