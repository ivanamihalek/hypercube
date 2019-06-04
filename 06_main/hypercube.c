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

Options options;
Keyword kwd[NO_KWDS];

int main ( int argc, char * argv[]) {

        Protein protein;
        Alignment alignment;
        int retval, method;
        int no_cons_methods, no_specs_methods;
        int *prot2almt = NULL, *almt2prot = NULL;
        int scores_assigned;
        int no_extern_spec_methods = 0;
        char ** score_label;
        double **spec_score  = NULL, **cons_score  = NULL;
        double *** in_group_score  = NULL;
        /********************************/
        int t, no_timesteps = 1000;
        double time_step    = 0.005;
        double ** rate_sym = NULL, *stnry_freq = NULL;
        double ** sim_matrix = NULL;
        double *** prob_matrix = NULL;

        /*******************************************/
        /*                                         */
        /*  input negotiation                      */
        /*                                         */
        /*******************************************/
        if (set_keywords ()) { /* (needed in the usage statement) */
                fprintf (stderr, "Error setting kwds (oink!?).\n");
                return 1;
        }

        /*****************/
        /* set defaults: */
        set_default_options ();
        set_aa2index();

        /* command file is required */
        if ( usage_statement (argc, argv)) return 1;

        /********************************************/
        /* change default options with the opt file */
        if  (read_opt_file (argv[1])) return 1;


        /*******************************************/
        /*                                         */
        /*  alignment -  input  and processing     */
        /*                                         */
        /*******************************************/

        /* read in the alignment */
        retval = read_alignment (&alignment);
        if (retval) exit(retval);


        /************************************************/
        /* do we know the protein groups in this story? */
        /* must happen before patching, because we are going
           to patch only within the group */
        if (options.groups[0]) {
                if (read_groups (options.groups, &alignment)) exit (1);
        } else {
                fprintf (stderr, "This implementation needs group definitions.\n");
                exit (1);
        }

        /*****************************************************************/
        /*   find  number of gaps, seq distances, "patch" the alignment  */
        retval =  process_almt( &alignment);
        if (retval) exit(retval);


        /*******************************************/
        /*                                         */
        /*  mapping on the structure, if provided  */
        /*                                         */
        /*******************************************/
        if ( options.pdbname[0]) {

                /* read in the structure */
                retval = read_pdb (options.pdbname, &protein, options.chain);
                if (retval) exit(retval);

                /* find mapping between the structure and the alignment*/
                if ( !(prot2almt = (int *) emalloc (protein.length*sizeof(int))) ) exit (1);
                if ( !(almt2prot = (int *) emalloc (alignment.length*sizeof(int))) ) exit (1);
                retval    = struct_almt_mapping (&protein, &alignment, prot2almt, almt2prot);
                if (retval) exit(retval);

        }
        /*******************************************/
        /*                                         */
        /*  surface, if provided                   */
        /*                                         */
        /*******************************************/
        if ( options.dssp_name[0] ) {
                /* surface = 0; */
                /* let's not couple the two - rather, the surface
                   could be made protected, but that's yet to be implemented*/
                retval  = read_dssp (&options, &protein);
                if (retval) exit(retval);
        }

        /*******************************************/
        /*                                         */
        /*  transition probabilities               */
        /*                                         */
        /*******************************************/
        if ( options.rate_matrix_file[0]) { // <<< the name of the file containing the
                // similairty (trans. prob, rather, matrix)
                if (!(rate_sym   = dmatrix(no_of_aa_types,no_of_aa_types))) return 1;
                if (!(stnry_freq = emalloc(no_of_aa_types*sizeof(double)))) return 1;
                if ( read_rate_matrix (options.rate_matrix_file, rate_sym, stnry_freq)) exit (1);

                /* allocate space for prob_matrices */
                if ( !(prob_matrix = emalloc (no_timesteps*sizeof(double**)) ) ) return 1;
                for (t=0; t<no_timesteps; t++) {
                        prob_matrix[t] = dmatrix (no_of_aa_types, no_of_aa_types);
                        if ( !prob_matrix[t]) return 1;

                }
                rate2prob (rate_sym, stnry_freq, time_step, no_timesteps, prob_matrix);

        } else { /* make it diagonal */
                int aa;
                sim_matrix = dmatrix (no_of_aa_types, no_of_aa_types);
                if ( !sim_matrix ) return 1;
                for (aa=0; aa<no_of_aa_types; aa++) sim_matrix[aa][aa] = 1.0;
        }


        /*******************************************/
        /*                                         */
        /* score the realiability of the alignment */
        /*                                         */
        /*******************************************/
        //will have to come back to this at some point ...
        //if ( !(reliability = emalloc(alignment.length*sizeof(double)))) exit (1);
        //retval = reliability_scoring (&alignment, reliability);
        //if (retval) exit(retval);


        /*******************************************/
        /*                                         */
        /*  conservation scoring                   */
        /*                                         */
        /*******************************************/
        /*********************************************************/
        /* the number of methods we are going to deal with here: */
        no_cons_methods = NUMBER_OF_CONS_METHODS;
        if ( !(cons_score = dmatrix( no_cons_methods, alignment.length))) exit (1);
        if ( !(in_group_score = emalloc(no_cons_methods*sizeof(double**)) )) exit (1);
        for (method=0; method < no_cons_methods; method++) {
                if ( !(in_group_score[method] = dmatrix(alignment.no_groups, alignment.length))) exit (1);
        }
        retval = conservation_scoring (&options, &alignment, NULL, cons_score, in_group_score);
        if (retval) exit(retval);

        /*******************************************/
        /*                                         */
        /*  specificity scoring                    */
        /*                                         */
        /*******************************************/
        if ( options.compare ) {
                /* options.extern_spec_methods contains both the filenames and the method
                   labels, so the actual number of methods is  options.extern_spec_methods[0][0]/2 */
                if ( options.extern_spec_methods ) {
                        no_extern_spec_methods = (int) options.extern_spec_methods[0][0]/2;
                } else {
                        no_extern_spec_methods = 0;
                }
                no_specs_methods = MAX_NO_SCORES + no_extern_spec_methods;
        } else {
                no_specs_methods = 1 + alignment.no_groups; /* discriminants + determinants for each group */
        }
        if ( !(spec_score    = dmatrix (no_specs_methods+no_extern_spec_methods, alignment.length))) exit (1);
        if ( !(score_label   = chmatrix(no_specs_methods+no_extern_spec_methods, SHORTSTRING))) exit (1);

        retval = spec_detection (&alignment, prob_matrix, no_timesteps, stnry_freq,
                                 spec_score, score_label, &scores_assigned);
        if (retval) exit(retval);


        for (method=0; method < no_extern_spec_methods; method++) {
                /* the first argument is the file name */
                retval =  read_extern_spec_method ( options.extern_spec_methods[2*method+1],
                                                    &alignment, spec_score[scores_assigned+method]);
                if ( retval) exit (retval);

        }



        /*******************************************/
        /*                                         */
        /*  output                                 */
        /*                                         */
        /*******************************************/
        score_output (&alignment, &protein, almt2prot, cons_score, in_group_score, spec_score);
        if (options.compare) {
                for (method=1; method <= no_extern_spec_methods; method++) {
                        sprintf (score_label[scores_assigned+method-1],"%s", options.extern_spec_methods[2*method]);
                }
                comparative_score_output (&alignment, spec_score, score_label, scores_assigned+no_extern_spec_methods);
        }

        /*******************************************/
        /*                                         */
        /*  cleanup                                */
        /*                                         */
        /*******************************************/
        almt_shutdown (&alignment);


        free_cmatrix(score_label);

        free_dmatrix(cons_score);
        free_dmatrix(spec_score);

        free (stnry_freq);

        for (t=0; t<no_timesteps; t++) free_dmatrix(prob_matrix[t]);
        free (prob_matrix);

        for (method=0; method < no_cons_methods; method++) {
                free_dmatrix (in_group_score[method]);
        }
        free (in_group_score);

        if ( options.pdbname[0]) {
                free (almt2prot);
                free (prot2almt);
        }

        if (options.rate_matrix_file[0]) {
        } else{
                free_dmatrix(sim_matrix);
        }
        return 0;
}

/***************************************************/
/***************************************************/
/***************************************************/
/***************************************************/
int usage_statement (int argc, char * argv[]) {

        int echo_kwds ();

        if (argc ==2 && !strcmp (argv[1], "-options" ) ) {
                /* if argc ==2 and the argv[1] is options
                   print out the available options */
                fprintf ( stderr, "Options, followed by the default value:\n" );
                echo_kwds();
                return 1;
        } else if (argc !=2 )   {
                /* otherwise print generic usage statement */
                fprintf ( stderr,
                          "Usage:  %s <options file>.\n",
                          argv[0]);
                fprintf ( stderr,
                          "   or:  %s -options \n\tto see the list of available options.\n",
                          argv[0]);
                return 1;
        }

        return 0;
}

/***************************************************/
/***************************************************/
int set_default_options () {

        memset (&options, 0, sizeof(Options) );
        /* what's this for? */

        options.max_gaps          =  0.3;
        options.patch_sim_cutoff  = -1.0;/* don't do the patching */
        options.patch_min_length  =  0.9;
        options.tree_method       = UPGMA;
        options.acc_cutoff        = 10.0;/* solvent acessibility cutoff */

        sprintf ( options.almtname, "%s", "in.msf");
        sprintf ( options.outname,  "%s", "cube");

        return 0;
}
/***************************************************/
/***************************************************/
int  grouping_max_strlen (Alignment *alignemnt) {

        int max;
        int nog = alignemnt->no_groups;
        int no_dec_places = 5;



        if ( nog >= 10000 ) {
                fprintf (stderr, "Well. %d might a few groups too many.\n", nog);
                exit (1);
        } else if (nog >= 1000 ) {
                no_dec_places = 4;
        } else if (nog >= 100  ) {
                no_dec_places = 3;
        } else if (nog >= 10   ) {
                no_dec_places = 2;
        } else {
                no_dec_places = 1;
        }

        max = (no_dec_places + 3)*nog; // 3 for ():

        return max;
}
