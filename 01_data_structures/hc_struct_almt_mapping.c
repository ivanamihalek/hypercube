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

int  struct_almt_mapping (Protein * protein, Alignment * alignment,
			  int * prot2almt, int * almt2prot){
    int prot_pos, almt_pos;
    char * struct_seq;
    Residue * prot_seq;
    /* locate query in the alignment */
    struct_seq = alignment->struct_seq;

    /*compare */
    prot_pos = 0;
    prot_seq = protein->sequence;
    for (almt_pos=0; almt_pos < alignment->length; almt_pos++ ) {
	
	if ( struct_seq [almt_pos] == '.' ) {
	    almt2prot [almt_pos] = -1;
	} else {
	    if ( prot_seq[prot_pos].res_type_short ==  struct_seq [almt_pos] ) {
		prot2almt[prot_pos] = almt_pos;
		almt2prot[almt_pos] = prot_pos;
	    } else {
		fprintf (stderr, "Structure/alignment mismatch,\n");
		fprintf (stderr, "\t structure: pdbid %s  value %c \n",
			 prot_seq[prot_pos].pdb_id,  prot_seq[prot_pos].res_type_short);
		fprintf (stderr, "\t alignment: pos %d  value  %c \n",   almt_pos,  struct_seq[almt_pos]);
		return 1;
	    }
	    prot_pos++;
	}
    }
    return 0;
}
