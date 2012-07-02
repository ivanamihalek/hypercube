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
