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
#ifndef PDB_H
# define	PDB_H



/* field start and length in atom record: */
    
/*  1 -  6        Record name     "ATOM  " */
# define        PDB_ATOM_REC_NAME          0
# define        PDB_ATOM_REC_NAME_LEN      6
/*     7 - 11        Integer         serial        Atom serial number. */
#define         PDB_ATOM_ATOM_NO           6       
#define         PDB_ATOM_ATOM_NO_LEN       5
  
/* 13 - 16        Atom            name          Atom name. */
#define         PDB_ATOM_ATOM_NAME         12
#define         PDB_ATOM_ATOM_NAME_LEN     4     
    

/* 17             Character       altLoc        Alternate location indicator. */
#define         PDB_ATOM_ALTLOC            16  
#define         PDB_ATOM_ALTLOC_LEN        1  

/* 18 - 20        Residue name    resName       Residue name. */
# define        PDB_ATOM_RES_NAME          17
# define        PDB_ATOM_RES_NAME_LEN      3

/* 22             Character       chainID       Chain identifier. */
#define        PDB_ATOM_CHAINID            21
#define        PDB_ATOM_CHAINID_LEN         1

/* 23 - 26        Integer         resSeq        Residue sequence number */
#define        PDB_ATOM_RES_NO             22
#define        PDB_ATOM_RES_NO_LEN         4

/* 27             AChar           iCode         Code for insertion of residues. */
#define        PDB_ATOM_INS_CODE           26
#define        PDB_ATOM_INS_CODE_LEN       1
    
/* 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in */
/*                                              Angstroms. */
#define        PDB_ATOM_X                  30
#define        PDB_ATOM_X_LEN               8
    
/* 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in */
/*                                              Angstroms. */
#define        PDB_ATOM_Y                  38
#define        PDB_ATOM_Y_LEN               8
    

/* 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in */
/*                                              Angstroms. */
#define        PDB_ATOM_Z                  46
#define        PDB_ATOM_Z_LEN               8
    


/* 55 - 60        Real(6.2)       occupancy     Occupancy. */
#define        PDB_ATOM_OCC                54
#define        PDB_ATOM_OCC_LEN             6
    

/* 61 - 66        Real(6.2)       tempFactor    Temperature factor. */
#define        PDB_ATOM_TEMP              60
#define        PDB_ATOM_TEMP_LEN           6  
    

/* 73 - 76        LString(4)      segID         Segment identifier, left-justified. */
#define        PDB_ATOM_SEGID             72
#define        PDB_ATOM_SEGID_LEN          4
    

/* 77 - 78        LString(2)      element       Element symbol, right-justified. */
#define        PDB_ATOM_ELMT              76
#define        PDB_ATOM_ELMT_LEN           2
     

/* 79 - 80        LString(2)      charge        Charge on the atom. */
#define        PDB_ATOM_CHARGE            78
#define        PDB_ATOM__CHARGE_LEN        2  
    

#endif /* PDB_H */
