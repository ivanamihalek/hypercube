# define ROOT  1
# define INNER 2
# define LEAF  4

# define UP    1
# define DOWN  2
# define LEFT  4
# define RIGHT 8

# define ASCII 128 /* number of characters */



typedef struct Node {

    struct Node *left, *right, *parent;
    int id;
    int type;
    int number_of_leaves;  /* in the corresponding subtree */
    int marked;
    int  bin; /*similarity bin the node belongs to*/
    double dist_to_parent, dist_to_left, dist_to_right;
    double avg_sim;
    char * seq;
    char * name;
    char consensus;
    int consensus_length;
    int population[ASCII]; /* population of amino acid types */
    double entropy;
    double entropy_complement;

} Node;

typedef struct Tree{

    Node *root;
    Node *leaf; /* this is actually the beginning of the node storage array */
    Node  ***group_root;   /* pointers to nodes which form groups at given rank */
    int size; /* leaves + inner nodes */
    int no_of_leaves;

} Tree;

