# define ALMT_NAME_LENGTH 50

typedef struct{
    char group_name[MEDSTRING];
    int  no_members;
    int  * group_member;
    int  *gaps; /* number of gaps at each position, in this group */
    int  *unks; /* number of unknowns at each position, in this group */
} Group;
    
typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int  *  seq_gaps;
    int  *  column_gaps;
    int  *  seq_unks;
    int  *  column_unks;
    int  * sunk; /* positions "sunk" to the bottom bcs of  too many gaps */
    int number_of_sunk_positions;
    double **seq_dist;
    int  **aligned_sites;
    int  **identical_sites;
    int  **similar_sites;
    char * struct_seq;
    int  * tp_flag; /* positions in the lmt tagged
		      for some reason as "true positive" */

    int no_groups;
    Group * group;
    int   * belongs_to_group;

    double *** trust_score;

    double ** per_node_spec_score;
    int    ** per_node_pool_size;
    
}  Alignment;
