# include "hypercube.h"

int  read_groups (char * filename, Alignment * alignment)  {
    FILE * fptr, *log = stdout;
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    char comment_char;
    int  line_ctr, retval;
    int  max_token;
    /***************/
    
    char **group_name, ***group_member;
    int group_ctr, member_ctr, i;
    int number_of_groups, no_of_members;
    int name_found;
    /***************/
    
    int echo_kwds ();
    int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
		 char * fmt, char * warnstr);
    /***************/
    
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;

    group_name = chmatrix (alignment->number_of_seqs, MEDSTRING);
    if ( !group_name) return 1;
    
    group_member = strmatrix (alignment->number_of_seqs,
			      alignment->number_of_seqs, ALMT_NAME_LENGTH);
    if ( !group_member) return 1;
    

    group_ctr = -1;
    member_ctr = -1;
    line_ctr = 0;
    memset ( line, 0, LONGSTRING);
    while(fgets(line,LONGSTRING,fptr)!=NULL){
 	line_ctr++;
 	/* tokenize */
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	switch ( retval ) {
	case  TOK_TOOMNY:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Too many tokens.");
	    fclose (log);
	    break;
	case TOK_TOOLONG:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Token too long.");
	    fclose (log);
	    break;
	}
	if ( max_token < 0 ) continue;
	
	
	/* if we are not reading a group we are looking for the kwd "name"*/
	if ( ! strcmp (token[0], "name")) {
	    
	    if ( max_token < 1 ) {
		errmsg ( log, line_ctr, line,
			 "\tKeyord %s should be followed by the group name.\n",
			 token[0]);
		return 1;
	    }
	    group_ctr  ++;
	    sprintf ( group_name[group_ctr], "%s", token[1]);
	    member_ctr = -1;
	} else if ( max_token == 0 ){

	    if ( group_ctr < 0 ) {
		errmsg ( log, line_ctr, line,
			 "\tEch group must be preceded by the keyword \"name\" followed by the group name.\n", "");
		return 1;
	    }
	    member_ctr ++;
	    sprintf ( group_member[group_ctr][member_ctr], "%s", token[0] );
	    
	} else if ( max_token == 0 ){
	    errmsg ( log, line_ctr, line,
		     "\tGroup members are expected one per line (change read_groups() if not happy).\n", "");
	    return 1;
	}
    }
    fclose (fptr);
    number_of_groups = group_ctr+1;

    /*********************************************************/
    /*allocate space for group info on the alignment pointer */
    alignment->no_groups = number_of_groups;
    alignment->group = emalloc (number_of_groups*sizeof(Group));
    /* for each group: */
    for (group_ctr = 0; group_ctr < number_of_groups; group_ctr++) {
	/* store group name */
	sprintf (alignment->group[group_ctr].group_name, "%s", group_name[group_ctr] );
	
	/* how many members in this group?*/
	member_ctr = -1;
	while ( (++member_ctr)<alignment->number_of_seqs
		&& group_member[group_ctr][member_ctr][0] );
	no_of_members = member_ctr;
	alignment->group[group_ctr].no_members = no_of_members;
	
	/* space for members */
	alignment->group[group_ctr].group_member = emalloc (no_of_members*sizeof(int) );
	if ( !alignment->group[group_ctr].group_member) return 1;

	
    }
    /* space for belongs-to array */
    alignment->belongs_to_group =  emalloc (alignment->number_of_seqs*sizeof(int) );
    if ( !alignment->belongs_to_group ) return 1;
    for (i=0; i<alignment->number_of_seqs; i++) alignment->belongs_to_group[i] = -1;

 
    /*********************************************************/
    /* sort out the group info                               */
    
    /* find the numerical id of each group member            */
    for (group_ctr = 0; group_ctr < number_of_groups; group_ctr++) {

	for ( member_ctr = 0; member_ctr< alignment->group[group_ctr].no_members; member_ctr++) {
	    name_found = 0;
	    for (i=0; i<alignment->number_of_seqs; i++) {
		if ( !strcmp (group_member[group_ctr][member_ctr], alignment->name[i]) ) {
		    alignment->group[group_ctr].group_member[member_ctr] = i;
		    alignment->belongs_to_group[i] = group_ctr;
		    name_found = 1;
		    break;
		}
	    }
	    if ( ! name_found ) {
		fprintf (stderr, "Name %s, listed in group %s in %s, not found in the alignment\n",
			 group_member[group_ctr][member_ctr], filename,
			 alignment->group[group_ctr].group_name);
		return 1;
	    }

	}
    }

    /* check if every sequence from the alignment belongs to some group */
    for (i=0; i<alignment->number_of_seqs; i++){
	if ( alignment->belongs_to_group[i] == -1 ) {
	    fprintf (stderr,
		     "Sequence %s from the alignment is not associated with any group given in %s.\n",
		     alignment->name[i], filename);
	    return 1;
	}
    }

    /* free the local space for groups */
    free_cmatrix   (group_name);
    free_strmatrix (group_member);

    return 0;
    

}





/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* junkyard  */

# if 0
    for (group_ctr = 0; group_ctr < number_of_groups; group_ctr++) {
	printf ("%d %s\n", group_ctr, alignment->group[group_ctr].group_name);
 	for ( member_ctr = 0; member_ctr< alignment->group[group_ctr].no_members; member_ctr++) {
	    i = alignment->group[group_ctr].group_member[member_ctr];
	    printf ("\t %s belongs to group %d\n", alignment->name[i], alignment->belongs_to_group[i]);
	}
    }
    exit (1);
# endif
# if 0
    printf ("**************************\n");
    for (i=0; i<alignment->number_of_seqs; i++){
	group_ctr = alignment->belongs_to_group[i];
	printf ("%s belongs to group %d (%s) \n", alignment->name[i],
		group_ctr, alignment->group[group_ctr].group_name);
    }
   
    //exit (1);
# endif
    
