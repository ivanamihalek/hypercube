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


typedef enum {OPTN_FLAG, OPTN_INTEGER, OPTN_DOUBLE, OPTN_EXP, OPTN_STRING, OPTN_INT_ARRAY, OPTN_STRING_ARRAY} Option_type;


/******************************/
/* tree building methods      */
/******************************/
# define UPGMA 1
# define NEIGHBOR_JOINING 2
# define CONSENSUS_UPGMA 3


/******************************/
/*   user options:            */
/******************************/
/* update the NO_KWDS above   */
/* and set_keywords () in cube_read_opt_file.c */

# define KWD_SIZE 30
# define NO_KWDS 20
# define MAX_STRING_ARRAY 10 /* for now used only for extern methods */

extern int no_kwds;

typedef struct {
    char pdbname [BUFFLEN];
    char struct_name [BUFFLEN]; /* the name (in the alignment)
				   of the sequence corresponding to the pdb structure */
    char dssp_name[BUFFLEN];  /* then name of the dssp file */
    char almtname[BUFFLEN];
    char outname [BUFFLEN];
    char groups  [BUFFLEN]; /* grouping of seqeunces into
			       presumptive functional or evolutionary groups */
    char rate_matrix_file[BUFFLEN]; /* file with the similarity matrix */
    char chain;
    int tree_score;
    double max_gaps;
    double patch_sim_cutoff;
    double patch_min_length;
    double acc_cutoff;  /* cutoff to call something solvent acessible */

    int guess_bad_translation;

    /* raw score vs. coverage: */
    int raw_score;
    int compare;  /* for papers and such - compare with other methods */
    char **extern_spec_methods; /*MAX_STRING_ARRAYxBUFFLEN char array will be allocated, if needed */

    int exchangeability; /* do something to account for the
			    exchangeability of the aa types (on by default if compare is on) */
    int tree_method;
    
} Options;

typedef struct {
    char name [KWD_SIZE];
    Option_type type;
    void * storage;
} Keyword;

/************************************/
extern Options options;       /* practically everybody needs acess to options */
extern Keyword kwd[NO_KWDS];  /* opt parser and main need to know about it */


    
