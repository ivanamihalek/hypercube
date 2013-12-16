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

int  read_dssp (Options * options, Protein *protein) {

    char line[LONGSTRING] = {'\0'};
    char pdb_id[PDB_ATOM_RES_NO_LEN+2] = {'\0'};
    char tmp[5];
    int resctr, acc, total = 0;
    FILE *fptr;
    
    fptr = efopen(options->dssp_name,"r");

    while( fgets( line, LONGSTRING, fptr) != NULL){
	if( ! strncmp(line+5,"RESIDUE", 7)){
	    break;
	}
    }
    while( fgets( line, LONGSTRING, fptr) != NULL){
	if ( strchr ( line, '!') ) continue;
	strncpy ( pdb_id, line+6, 5);
	string_clean (pdb_id, 5);
	if ( ! pdb_id[0] )  continue;
	strncpy ( tmp, line+34, 4);
	tmp[4] = '\0';
	acc = atoi(tmp);
	if ( acc > options->acc_cutoff ) {
	    for ( resctr=0; resctr < protein->length; resctr++ ) {
		if ( !strcmp (pdb_id, protein->sequence[resctr].pdb_id) ) {
		    protein->sequence[resctr].solvent_accessible = 1;
		    total++;
		}
	    }
	}
    }
    fclose (fptr);

    if ( ! total ) {
	fprintf (stderr, "No residue solvent accessible (?!)\n");
	return 1;
    }

    return 0;
}
