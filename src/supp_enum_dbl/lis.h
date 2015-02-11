#include <gmp.h>
#include "../utils/matrix.h"

#ifndef LIS_H
#define LIS_H

#define LIS_RREF 1
#define LIS_LU 2

int lis_rref_has_infinite(matrix_t *A);
int lis_rref_has_no_sol(matrix_t *A);
matrix_t* lis_solve(matrix_t *A, matrix_t *b);
matrix_t* lis_solve_system(matrix_t *A, matrix_t *b, int method, int *err);

#endif
