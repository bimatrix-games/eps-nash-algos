#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"

#ifndef IO_H
#define IO_H
#define MAXSTR 100

void not_impl (char *info);
void read_conf (const char *s);
matrix_t **read_bimatrix_from_file(FILE *f, int* rdim1, int* rdim2) ;
matrix_t **read_bimatrix();
matrix_t **create_matrices(int m, int n);
void free_matrices(matrix_t **mats);
#endif
