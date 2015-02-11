#include <math.h>
#include "../utils/matrix.h"

#ifndef UTIL_H
#define UTIL_H

typedef struct {
    matrix_t *xt;
    matrix_t *Ct;
    matrix_t *xt_R_y;
    matrix_t *xt_C_y;
    matrix_t *Ct_x;
    matrix_t *xt_C;
    matrix_t *R_y;
    int br_R_y_index;
    double br_R_y_val;
    int br_Ct_x_index;
    double br_Ct_x_val;
}bbm_extra_t;

bbm_extra_t *bbm_extra_alloc(matrix_t *R, matrix_t *C, matrix_t *x, matrix_t *y);
void bbm_extra_free(bbm_extra_t *ext);

double best_response(matrix_t *mat_s, int *index);
double compute_epsilon(bbm_extra_t *ext, double *eps1, double *eps2);
#endif
