#ifndef SYMNMF_H
#define SYMNMF_H

#include <stdlib.h>

/* Function prototypes */
double** sym(double **data, int N, int d);
double** ddg(int N, double **A);
double** norm (int N, double **D, double **A);
void symnmf(double **curr_H, double **W, int N, int d); 

void allocate_matrix_memory(double ***matrix, int N, int k);
void free_matrix_memory(double **matrix, int rows);

#endif /* SYMNMF_H */
