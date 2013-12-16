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
# ifndef _UTILS_H
# define _UTILS_H
# include <stdio.h>

# define ASCII 128 

extern char *amino_acid_order;
extern int no_of_aa_types;
extern int aa2index[ASCII];


int      array_qsort (int * sorted_pos, double * sa, int sequence_length );
char   **chmatrix(int rows, int columns);
double **dmatrix(int rows, int columns);
double ***d3matrix(int rows, int columns, int floors);
char   ***strmatrix(int rows, int columns, int floors);
void *   emalloc(int	size);
FILE *   efopen(char * name, char * mode);
void     free_cmatrix(char **m);
void     free_imatrix(int **m);
void     free_dmatrix(double **m);
void     free_d3matrix(double ***m);
void     free_strmatrix(char ***m);
int    **intmatrix(int rows, int columns);
int set_aa2index ();
int      string_clean ( char* string, int length);
# endif
