# include "hypercube.h"

/********************************************************************************/
int build_tree (Options * options, Alignment * alignment, Tree * tree) {
    
# define NUMBER_OF_BINS 20

    int rank, group_tag, retval;
    int upgma_and_nj (Options * options, Alignment * alignment, Tree * tree);
    int find_group_roots  (Node *node, Node*** group_root,  int rank, int * group_tag );
    Node  ***node_matrix(int rows, int columns);
    
    switch (options->tree_method) {
    case UPGMA:
    case CONSENSUS_UPGMA:
    case NEIGHBOR_JOINING:
	retval =  upgma_and_nj (options, alignment, tree);
	 /* should fill the tree->size value*/
	if ( retval ) return retval;
	break;
    default:
	fprintf ( stderr,"Unrecognized tree building method.\n");
	return 1;
    }
    
    
    /* find groups */
    tree->group_root = node_matrix (tree->no_of_leaves, tree->no_of_leaves);
    for ( rank = 1; rank < tree->no_of_leaves; rank++ ) {
	group_tag = -1;
	find_group_roots (tree->root,  tree->group_root,  rank, & group_tag );
    }

    return 0;
}


/**********************************************************************************************/
int place_root ( Node *leaf, int no_of_nodes, Node * root) {
    
    int edge_ctr = 0, longest_ctr = 0, i;
    int find_path ( Node ** path, int * edge_ctr, Node ** longest_path,  int * longest_ctr,Node * node, Node* from) ;

    Node** path, ** longest_path;
    if ( ! (path = emalloc (no_of_nodes*sizeof(Node*) ) ) ) exit (1);
    if ( ! (longest_path = emalloc (no_of_nodes*sizeof(Node*) ) ) ) exit (1);
    /* print all possible paths: */
    /* find any leaf to start fom: */
    longest_ctr = 0;
    edge_ctr = 0;
    for (i = 0; i < no_of_nodes; i++ ) {
	if ( (leaf + i) -> type != LEAF  ) continue;
	find_path (path, &edge_ctr, longest_path, &longest_ctr, leaf + i, NULL);
   }
    /* printf (" %4d  %s --->  %s \n", longest_ctr, longest_path[0]->name,  longest_path[longest_ctr-1]->name); */
    
    /* place the root above the middle node on that path */
    {
	int insert_node ( Node * new_node, Node * old_node_1, Node * old_node_2);
	int reorient_nodes ( Node * node, Node * from );
	root->type = ROOT;
	insert_node (root, longest_path[longest_ctr/2-1], longest_path[longest_ctr/2]);
 	reorient_nodes ( root, NULL);
    }
    return 0;
}


/**********************************************************************************************/
int rank_inner_nodes (Node *root, int  no_of_nodes) {

    int i, inner_node_ctr;
    int  * node_label_sorted;
    double *distance;
    Node ** label2node;
    int find_leaves_dist ( Node *node, int *inner_node_ctr, int * node_label_sorted, Node ** label2node, double * distance );
    
    /* allocate */
    if ( ! ( node_label_sorted = emalloc (no_of_nodes*sizeof(int)))) return 1;
    if ( ! ( distance  = emalloc (no_of_nodes*sizeof(double)))) return 1;
    if ( ! ( label2node  = emalloc (no_of_nodes*sizeof(Node*)))) return 1;
    
    /*find distance to the root from each inner node */
    inner_node_ctr = 0;
    find_leaves_dist ( root, &inner_node_ctr, node_label_sorted, label2node, distance );
   

    /* sort inner nodes by that distance */
    array_qsort ( node_label_sorted, distance, inner_node_ctr);
    /* assign them their sorted number */
    for (i = 0; i < inner_node_ctr; i++ ) {
	label2node[node_label_sorted[i]]->id = i+1;
    }
   
    free (node_label_sorted);
    free (label2node);
    free (distance);

    return 0;
}
/********************************************************************************/
double  node_distance ( double **seq_dist, Node* node1, Node* node2 ){

    
    if (  (node1->type == LEAF) &&  (node2->type == LEAF) ) {
	return seq_dist[node1->id][node2->id];
    } else {
	double distance = 0.0;
	if ( node1->type != LEAF ) { /* left side depth first */
	    distance += node_distance( seq_dist, node1->left,  node2);
	    distance += node_distance( seq_dist, node1->right, node2);
	} else 	if ( node2->type != LEAF ) {
	    distance += node_distance( seq_dist, node1,  node2->left);
	    distance += node_distance( seq_dist, node1,  node2->right);
	}
	return distance;
    }
}
/**********************************************************************************************/
int insert_node ( Node * new_node, Node * old_node_1, Node * old_node_2) {
    new_node->left = old_node_1;
    new_node->right= old_node_2;

    if  ( old_node_1 ==  old_node_2->left ) {
	old_node_2->left = new_node;
    } else if  ( old_node_1 ==  old_node_2->right ){
	old_node_2->right = new_node;
    } else if  ( old_node_1 ==  old_node_2->parent ){
	old_node_2->parent = new_node;
	old_node_2->dist_to_parent = old_node_1->dist_to_parent =  old_node_2->dist_to_parent/2;
	
    } else {
  	fprintf ( stderr, "Error inserting node.\n");
	exit (1);
    }

    if  ( old_node_2 ==  old_node_1->left ) {
	old_node_1->left = new_node;
    } else if  ( old_node_2 ==  old_node_1->right ){
	old_node_1->right = new_node;
    } else if  ( old_node_2 ==  old_node_1->parent ){
	old_node_1->parent = new_node;
 	old_node_2->dist_to_parent = old_node_1->dist_to_parent =  old_node_1->dist_to_parent/2;
   } else {
  	fprintf ( stderr, "Error inserting node.\n");
	exit (1);
    }

    new_node->dist_to_left  =  old_node_1 -> dist_to_parent;
    new_node->dist_to_right =  old_node_2 -> dist_to_parent;
    return 0;
}


/**********************************************************************************************/
int reorient_nodes ( Node * node, Node * from ) {

    Node *left, *right;
    double dist_to_left, dist_to_right, dist_to_parent;
    /* decide on "left" and "right" */
    if (  from ==  node->parent  ) {
	
	left = node->left;
	right = node->right;
	
	dist_to_left   = node->dist_to_left;
	dist_to_right  = node->dist_to_right;
	dist_to_parent =  node->dist_to_parent;
	
    } else if (from  ==  node->left  ) {
	
	right = node->right;
	left  = node->parent;
	
	dist_to_right  = node->dist_to_right;
	dist_to_left   = node->dist_to_parent;
	dist_to_parent = node->dist_to_left;
	
    } else if ( from == node->right  ) {
	
	left   = node->left;
	right  = node->parent;
	
	dist_to_left   = node->dist_to_left;
	dist_to_right  = node->dist_to_parent;
	dist_to_parent = node->dist_to_right;

	
	
    } else {
	fprintf ( stderr, "Error traversing tree in reorient_nodes.\n"); 
	exit (1);
    }
    node->parent = from; 
    node->left   = left; 
    node->right  = right; 
    node->dist_to_left   = dist_to_left;
    node->dist_to_right  = dist_to_right;
    node->dist_to_parent = dist_to_parent;
 
    if ( node->type != LEAF ) { 
	reorient_nodes (left,  node); 
	reorient_nodes (right, node);
    }
    return 0;
}
/**********************************************************************************************/
int find_path ( Node ** path, int * edge_ctr, Node ** longest_path,  int * longest_ctr,Node * node, Node* from) {

    if ( ! node )  {
	fprintf ( stderr, "Error traversing tree in find_path (incomming node = 0x0).\n");
	exit (1);
    }
    
    /* add yourself to the path */
    path[*edge_ctr] = node;
    /* increase the counter     */
    (*edge_ctr) ++;
    if (  node->type != LEAF ) {
	Node *left, *right;
	/* decide on "left" and "right" */
	if (  node->parent == from ) {
	    left = node->left;
	    right = node->right;
	} else if ( node->left == from ) {
	    right = node->right;
	    left  = node->parent;
	} else if ( node->right == from ) {
	    right = node->parent;
	    left  = node->left;
	} else {
	    fprintf ( stderr, "Error traversing tree in find_path.\n");
	    exit (1);
	}
	/* go down the left */
	find_path (path, edge_ctr, longest_path, longest_ctr,  left, node);
	/* go down the   right */
	find_path (path, edge_ctr, longest_path, longest_ctr,  right, node);
    } else {
	if ( *edge_ctr == 1 ) { /* we are at the beginning */
	    find_path (path, edge_ctr, longest_path, longest_ctr,node->parent, node);
	} else {  /* we are at the end */
	    if ( *edge_ctr > *longest_ctr ) {
		*longest_ctr = *edge_ctr;
		memcpy ( longest_path, path, (*longest_ctr)*sizeof(Node*));
	    }
	}
    }
    /* get yourself off the path */
    path[*edge_ctr] = NULL;
    /* decrease the counter */
    (*edge_ctr) --;
  
    return 0;
}
/**********************************************************************************************/
int find_leaves_dist ( Node *node, int *inner_node_ctr, int * node_label_sorted, Node ** label2node, double * distance ){
    double distance_to_parent (Node * node);
 
    if ( node -> type == LEAF ) {
    } else {
	node_label_sorted[*inner_node_ctr] = *inner_node_ctr;
	label2node[*inner_node_ctr] = node;
	distance[*inner_node_ctr] = distance_to_parent (node);
 	(*inner_node_ctr) ++;
	find_leaves_dist ( node->left, inner_node_ctr, node_label_sorted, label2node, distance );
	find_leaves_dist ( node->right, inner_node_ctr, node_label_sorted, label2node, distance );
    }

    return 0;
}

/**********************************************************************************************/
double distance_to_parent (Node * node){

    Node * current = node ;
    double distance= 0.0;
    
    while ( current->parent ) {
	distance += current -> dist_to_parent;
	current = current->parent;
    };

    return distance;
    
}




/********************************************************************************/

/*******************************************************************/
Node  ***node_matrix(int rows, int columns){
    Node ***m;
    int i;
        /* allocate pointers to rows */
    m=(Node ***) malloc(rows*sizeof(Node**));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(Node **) calloc( rows*columns, sizeof(Node*));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

int free_node_matrix ( Node  *** m ) {
    free (m[0]);
    free (m);
    return 0;
}


/******************************************************************************/
int find_group_roots (Node *node, Node*** group_root,  int rank, int * group_tag ) {
    if ( node->type==LEAF ) {
    } else {
	
	/* determining the roots for all the groups at the rank "rank" */
	/* remember:  at rank N, N is the smallest group root possible */
	if ( node->id == rank ) {
	    ++(*group_tag);
	    group_root[rank][*group_tag] = node;
	} else if( node->id < rank ){
	    
	    Node *left, *right;
	    left  = node->left;
	    right = node->right;

	    if ( left->id >= rank  || left->type==LEAF) {
		++(*group_tag);
		group_root[rank][*group_tag] = left;
	    } else {
		find_group_roots ( left, group_root, rank, group_tag);
	    }
	    if ( right->id >= rank  || right->type==LEAF) {
		++(*group_tag);
		group_root[rank][*group_tag] = right;
	    } else {
		find_group_roots ( right, group_root, rank, group_tag);
	    }
	}
	
    }
    return 0;
}

/********************************************************************************/
/********************************************************************************/
int upgma_and_nj (Options * options, Alignment * alignment, Tree * tree) {

    double  distance1, distance2;
    int retval;
    int node_ctr, no_seqs = alignment->number_of_seqs, tree_size;
    int closest1, closest2,* current_list;
    int upper, last_ctr;
    Node * node;
    int ( *closest_in_curr_list ) (Alignment * alignment, Node * node, int tree_size, int new_node, 
				   int *  current_list, int * closest1, int *closest2,
				   double * dist_1_ptr, double * dist_2_ptr );
    int closest_in_curr_list_nj (Alignment *alignment, Node * node, int tree_size, int new_node, 
			     int *  current_list, int * closest1, int *closest2,
				 double * dist_1_ptr, double * dist_2_ptr );
   /* I need dummy here to make upgma and nj defs look formally the same */ 
    int closest_in_curr_list_upgma (Alignment *alignment, Node * node, int tree_size, int dummy, 
			     int *  current_list, int * closest1, int *closest2,
				    double * dist_1_ptr, double * dist_2_ptr );
    int fill_consensus_array (int seq_length, char * seq1, char * seq2, char *consensus);
    switch (options->tree_method) {
    case UPGMA:
	closest_in_curr_list = closest_in_curr_list_upgma;
	break;
   case NEIGHBOR_JOINING:
	closest_in_curr_list = closest_in_curr_list_nj;
	break;
    default:
	fprintf ( stderr,"Unrecognized tree building method.\n");
	return 1;
    }

    /*****************************************/
    /*****************************************/
    /*          build tree                   */
    /*****************************************/
    /*****************************************/
    
    /* allocate space */
    tree_size = 2*no_seqs -1;
    node = (Node *) emalloc ( tree_size*sizeof(Node) );
    /* initialize leaves to point to sequences from the alignment*/
    for ( node_ctr=0; node_ctr < no_seqs; node_ctr++ ) {
	node[node_ctr].id   = node_ctr;
	node[node_ctr].type = LEAF;
	node[node_ctr].seq  = alignment->sequence[node_ctr];
        node[node_ctr].name = alignment->name[node_ctr];
	node[node_ctr].number_of_leaves = 1;
    }
    /* initialize current list of nodes whose distance needs to be compared */
    current_list = (int*) emalloc ( tree_size*sizeof(int));
    for ( node_ctr=0; node_ctr < no_seqs; node_ctr++ ) {
	current_list[node_ctr] = 1; 
    }

    /* find children for each of the remaining nodes */
    upper = (options->tree_method == NEIGHBOR_JOINING) ? tree_size - 1 : tree_size;
    for ( node_ctr=no_seqs; node_ctr < upper; node_ctr++ ) {
	retval = closest_in_curr_list (alignment, node, tree_size, node_ctr,
				       current_list, &closest1, &closest2,
				       &distance1, &distance2);
	if (retval) return retval;
	/* hack, so I can order the nodes in the tree: */
	if ( distance1 <= 0 ) distance1 = 0.001;
	if ( distance2 <= 0 ) distance2 = 0.001;
	
	/* fill in the new  node fields */ 
	node[node_ctr].left             = &node[closest1];
	node[node_ctr].dist_to_left     = distance1;
	node[node_ctr].right            = &node[closest2];
	node[node_ctr].dist_to_right    = distance2;
	node[node_ctr].type             = INNER;
	node[node_ctr].id               = tree_size - node_ctr; /* this will serve as the "rank" */
	node[node_ctr].number_of_leaves = node[closest1].number_of_leaves +
	    node[closest2].number_of_leaves;
	if ( ! 	node[node_ctr].seq ) {
	    /* what do I need this for ? */
	    //node[node_ctr].seq   = emalloc ( alignment->length*sizeof(char) );
	    //if ( ! node[node_ctr].seq ) exit (1);
	}
	    
	node[closest1].parent = & node[node_ctr];
	node[closest2].parent = & node[node_ctr];
	node[closest1].dist_to_parent = distance1;
	node[closest2].dist_to_parent = distance2;
	/* remove the two from the current list, and replace them with the parent */ 
	current_list[closest1] = 0;
	current_list[closest2] = 0;
	current_list[node_ctr] = 1;
    }

    last_ctr = node_ctr - 1;
    tree->leaf = node;
    (node + tree_size -1)->type = ROOT;
    tree->root = node + tree_size - 1;

    
    if( options->tree_method == NEIGHBOR_JOINING) {
	
	int place_root ( Node *leaf, int no_of_nodes, Node * root) ;
	int rank_inner_nodes (Node *root, int  no_of_nodes);
	/* make the last two nodes  parents of each other */
	/* place the root */
	for ( node_ctr=0; node_ctr < upper; node_ctr++ ) {
	    if ( current_list[node_ctr] ) printf ( " %d ", node_ctr);
	}
	printf ("\n");
	for ( node_ctr=0; node_ctr < upper; node_ctr++ ) {
	    if ( current_list[node_ctr] ) {
		node[node_ctr].parent =  & node[last_ctr];
		node[last_ctr].parent =  & node[node_ctr];
		break;
	    }
	}
	place_root ( tree->leaf,  tree_size-1, tree->root);
	rank_inner_nodes ( tree->root,  tree_size);
    }

    tree->size = tree_size;
    tree->no_of_leaves = no_seqs; 
# if 0
    print_tree (stdout, tree->root);
    exit(0);
# endif
    
    free (current_list);
   
    return 0;

}
/********************************************************************************/
/********************************************************************************/
int closest_in_curr_list_upgma (Alignment * alignment, Node * node, int tree_size, int dummy, 
				int *  current_list, int * closest1, int *closest2,
				double * dist_1_ptr, double * dist_2_ptr ){
    int ctr1, ctr2;
    double distance, min_distance;
    double **seq_dist = alignment->seq_dist;
    double  node_distance ( double **seq_dist, Node* node1, Node* node2 );

    min_distance = 1000;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
        if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    distance = node_distance ( seq_dist,  node+ctr1, node+ctr2 );
	    //exit (0);
	    /* we need the average */ 
	    distance /= (node+ctr1)->number_of_leaves*(node+ctr2)->number_of_leaves;
	    if ( distance > 1 ) {
		printf ( "closest_in_curr_list_upgma(): %d   %d  %d   %d  %8.2lf\n",
			 ctr1, ctr2,  (node+ctr1)->number_of_leaves,
			 (node+ctr2)->number_of_leaves, distance);
		almt_shutdown (alignment);
		exit (1);
	    }
	    if ( distance < min_distance) {
		min_distance = distance;
		*closest1 = ctr1;
		*closest2 = ctr2;
	    }
	}
    }
    *dist_1_ptr = min_distance/2;
    *dist_2_ptr = min_distance/2;
    return 0;
}
int closest_in_curr_list_nj (Alignment *alignment, Node * node, int tree_size, int new_node, 
			     int *  current_list, int * closest1, int *closest2,
			     double * dist_1_ptr, double * dist_2_ptr ){
    /* actually, for nj, they are not the closest, but rather the pair which,
       if joined next, gives the minimum total sum of branch lengths */
    /* new_node  is where I'll put the parent for the closest pair I find */
    
    int ctr1, ctr2, curr_list_length;
    int no_leaves = (tree_size+1)/2;
    double sum, min_sum;
    double avg1, avg2, avg1_min, avg2_min;
    double **seq_dist = alignment->seq_dist;
    static double ** dist_table = NULL;
    double  nj_sum_of_branch_lengths (double ** dist_table,  int * current_list,
				   int  node_ctr1, int node_ctr2,
				   int tree_size, double *avg1_ptr, double *avg2_ptr);

    if ( ! dist_table ) { /* initialize distance table */
	/* allocate */
	dist_table = dmatrix (tree_size, tree_size);
	/* copy distances btw the leaves */
	for (ctr1=0; ctr1 < no_leaves; ctr1++ ) {
	    for (ctr2=ctr1+1; ctr2 < no_leaves; ctr2++ ) {
		dist_table[ctr1][ctr2] = dist_table[ctr2][ctr1] = seq_dist[ctr1][ctr2];
	    }
	}
    }
    
    curr_list_length = 0;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	curr_list_length += current_list[ctr1];
    }

    avg1_min = 0;
    avg2_min = 0;
    min_sum = 10000;
    
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    sum = nj_sum_of_branch_lengths ( dist_table, current_list,
					     ctr1, ctr2, tree_size, &avg1, &avg2);
	    if ( sum < min_sum) {
		min_sum   = sum;
		*closest1 = ctr1;
		*closest2 = ctr2;
		avg1_min  = avg1; 
		avg2_min  = avg2; 
	    }
	}
    }
    
    printf ( "closest:  %3d  %3d  \n\n", *closest1, *closest2);
    
    *dist_1_ptr = ( dist_table[*closest1][*closest2]+ avg1_min - avg2_min)/2;
    *dist_2_ptr = ( dist_table[*closest1][*closest2]+ avg2_min - avg1_min)/2;

    dist_table[new_node][*closest1] =  dist_table[*closest1][new_node] = *dist_1_ptr;
    dist_table[new_node][*closest2] =  dist_table[*closest2][new_node] = *dist_2_ptr;

    /* update the distance table */
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	if ( ctr1 == new_node )   continue;
	if ( ctr1 == *closest1 ||  ctr1 == *closest2) continue;
	dist_table[new_node][ctr1] = dist_table[ctr1][new_node] =
	    ( dist_table[*closest1][ctr1] + dist_table[*closest2][ctr1]
	      -  dist_table[*closest2][*closest1] )* 0.5;
    }
	
	
 
    if ( curr_list_length == 3 ) {
	for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	    for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
		printf ( " %3d  %3d  %8.4lf \n", ctr1, ctr2,  dist_table[ctr1][ctr2]);
	    }
	    printf ( "\n");
	}
    }
    /* printf (" %4d    %4d  %4d   %8.3lf   %8.3lf   %8.3lf   %8.3lf   %8.3lf  %8.3lf   \n", */
    /*  curr_list_length,  *closest1, *closest2, dist_table[*closest1][*closest2],
	min_distance, avg1_min, avg2_min, *dist_1_ptr, *dist_2_ptr); */
   
    return 0;
}


/***************************************************************************************/
double  nj_sum_of_branch_lengths (double ** dist_table,  int * current_list,
				  int  node_ctr1, int node_ctr2, int tree_size,
				  double *avg1_ptr, double *avg2_ptr) {
    double sum, term;
    double avg[2];
    int current_pair[2] = {node_ctr1,node_ctr2};
    int ctr1, ctr2, n;
    
    sum =  dist_table[node_ctr1][node_ctr2]/2;


    for (ctr1 = 0; ctr1 < 2; ctr1++ ) {
	n =0;
	avg[ctr1] = 0;
	for (ctr2= 0; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    if ( ctr2 == node_ctr1 || ctr2  == node_ctr2 )
		continue;
	    avg[ctr1] += dist_table[ current_pair[ctr1]][ ctr2];
	    n ++;
	}
	if ( !n ) {
	    printf ( "** %d %d \n", node_ctr1, node_ctr2);
	    printf ( "available:  \n");
	    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
		if ( ! current_list[ctr1]) continue;
		printf ( "%4d", ctr1);
	    }
	    printf ( "\n");
	    exit (1);
	}
    
 	    
	avg[ctr1] /= n;
    }

    sum += (avg[0] + avg[1])/2;
    
    *avg1_ptr = avg[0];
    *avg2_ptr = avg[1];
    
    term = 0;
    n = 0;
    for (ctr1=0; ctr1 < tree_size; ctr1++ ) {
	if ( ! current_list[ctr1]) continue;
	if ( ctr1 == node_ctr1 ||  ctr2  == node_ctr2 )
	    continue;
	for (ctr2=ctr1+1; ctr2 < tree_size; ctr2++ ) {
	    if ( ! current_list[ctr2]) continue;
	    if ( ctr2 == node_ctr1 ||  ctr2 == node_ctr2 )
		continue;
	    term += dist_table[ctr1][ctr2];
	    n ++;
	}
    }
    sum += term/n;
		
	    
    return sum;
}
