/*******************************************/
/* by  Ivana Mihalek, 2011                 */
/*******************************************/

# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>

# define BUFFLEN     150
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25

/******************************/
/*   tokenizer                */
/******************************/

# define TOK_TOOMNY  1 /* tokenizer error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */

# include "hc_pdb.h"
# include "hc_alignment.h"
# include "hc_geometry.h"
# include "hc_utils.h"
# include "hc_options.h"
# include "hc_tree.h"

# define VERSION "2009"


/*****************************/
/*   log   modes             */
/*****************************/
# define INTRO 1
# define WARN  2
# define STATUS  3
# define NOTE  4


/******************************/
/*   geometry                 */
/******************************/
/* distance to be called a neighbor in the protein structure */
# define CUTOFF_DIST 4.0




/******************************/
/*  similarity scores         */
/******************************/
# define NUMBER_OF_CONS_METHODS  4  /* entropy and rvet for now */


typedef enum {DISCR, DET, NUMBER_OF_SPECS_METHODS} Method_type;
extern char method_name[NUMBER_OF_SPECS_METHODS][SHORTSTRING];

# define MAX_NO_SCORES 100

# define ALL_GAPS  1000
# define HALF_WIN 2  /* used to estimate the alignment reliabilty */




/******************************/
/* function declarations :    */
/******************************/

int afa_out (Options *options, Alignment * alignment);
int almt_shutdown  (Alignment *alignment);
int alt_spec_detection (Alignment * alignment, double ** similarity, double * stnry_freq, double **score);

void cluster_counter (int  no_of_things,  int *neighbors[],  int * clusters[]); 
int coverage ( Alignment * alignment, double * score,  double * res_fract_rank, int group_id);

int tokenize ( char  token[MAX_TOK][MEDSTRING], int * max_token,
	       char * line, char comment_char);

int chisquare ( int * population_1, int * population_2, int no_bins, int kstr,
		double *df, double *chsq, double *prob);
int comparative_score_output (Alignment *alignment, double ** score, char**score_label, int no_scores);
int conservation_scoring (Options *options, Alignment * alignment, 
			  int * similar_to, double ** score, double *** in_group_score);
int cube_output   (Alignment * alignment, double  **rate_sym, double * freq, double * score);
int free_node_matrix (Node  ***m);
int normalize (double *array, int size, int power);
int patch_almt    (Alignment * alignment);
int process_almt  (Alignment * alignment);
     
int rate2prob (double ** rate_sym, double * freq, double time_step, int no_timesteps, double *** prob_matrix);
int read_clustalw (Alignment * alignment);
int read_extern_spec_method ( char *filename,  Alignment * alignment, double *spec_score);
int read_groups (char * filename, Alignment * alignment);
int read_opt_file (char * filename );
int read_pdb ( char * pdbname, Protein * protein, char chain);
int read_rate_matrix (char * infile_name, double **rate_sym, double *freq) ;

int rocs_output (Alignment *alignment, double ** score);
int score_almt_trust (Alignment *alignment, int similar_to[]);
int score_output (Alignment *alignment, Protein *protein, int *almt2prot, 
		  double ** score, double *** in_group_score, double ** spec_score);  
int score_trees (Alignment * alignment);
int scoring ( Alignment * alignment, double sim_matrix[][ASCII],
	      double ** score, double *reliability);
int seq_similarity_indicators (Alignment * alignment);
int set_default_options ();
int set_keywords ();

int spec_detection ( Alignment * alignment,  double *** prob_matrix_series, int no_timesteps,
		      double * stnry_freq, double **spec_score, char ** score_label, int *scores_assigned);
int struct_almt_mapping (Protein * protein, Alignment * alignment,
			 int * prot2almt, int * almt2prot);
int usage_statement (int argc, char * argv[]);
 



# endif
