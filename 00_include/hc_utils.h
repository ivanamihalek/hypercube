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
