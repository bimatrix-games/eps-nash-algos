#include "util.h"

/* Allocates and computes extra resources neccessary for computation */
bbm_extra_t *bbm_extra_alloc(matrix_t *R, matrix_t *C, matrix_t *x, matrix_t *y)
{
    bbm_extra_t *ext = malloc(sizeof(bbm_extra_t));

    ext->xt = matrix_trans(x);
    ext->Ct = matrix_trans(C);
    ext->Ct_x = matrix_mul_mat_vec(ext->Ct, x);
    ext->xt_C = matrix_trans(ext->Ct_x);
    ext->R_y = matrix_mul_mat_vec(R, y);
    ext->xt_R_y = matrix_mul_vec_mat(ext->xt, ext->R_y);
    ext->xt_C_y = matrix_mul_mat_vec(ext->xt_C, y);
    ext->br_R_y_val = best_response(ext->R_y, &(ext->br_R_y_index));
    ext->br_Ct_x_val = best_response(ext->Ct_x, &(ext->br_Ct_x_index));

    return ext;
}

/* Frees up allocated resources */
void bbm_extra_free(bbm_extra_t *ext)
{
    matrix_free(ext->xt);
    matrix_free(ext->Ct);
    matrix_free(ext->xt_R_y);
    matrix_free(ext->xt_C_y);
    matrix_free(ext->Ct_x);
    matrix_free(ext->xt_C);
    matrix_free(ext->R_y);
    free(ext);
}

/* 
 * Computes the best response given a vector of expected payoffs
 * for each pure strategy.
 */
double best_response(matrix_t *mat_s, int *index)
{
    int i, j;
    double payoff = mat_s->data[0][0];

    j = 0;
    for (i = 0; i < mat_s->nrows; ++i)
        if (payoff < mat_s->data[i][0]) {
            j = i;
            payoff = mat_s->data[i][0];
        }

    if (index != NULL)
        *index = j;

    return payoff;
}

/* 
 * Computes the epsilon-NE for both players and stores them in eps1 and eps2
 * for player 1 and 2 respectively and returns the maximum of the two.
 */
double compute_epsilon(bbm_extra_t *ext, double *eps1, double *eps2)
{
    double u, v, u1, v1, e1, e2;
    
    u = ext->xt_R_y->data[0][0];
    v = ext->xt_C_y->data[0][0];
    u1 = ext->br_R_y_val;
    v1 = ext->br_Ct_x_val;
    
    e1 = (u < u1) ? u1 - u : 0;
    e2 = (v < v1) ? v1 - v : 0;

    if (eps1 != NULL)
        *eps1 = e1;
    if (eps2 != NULL)
        *eps2 = e2;

    if (e1 > e2)
        return e1;

    return e2;
}
