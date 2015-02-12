#include <math.h>
#include "checker.h"
#include "lis.h"
#include "../utils/util.h"

matrix_t *R;
matrix_t *C;
matrix_t *x;
matrix_t *y;
matrix_t *supp_x;
matrix_t *supp_y;

/* Form the linear system given a players payoff matrix */
matrix_t** form_lis(matrix_t *A)
{
    int i, j, k, m, n;
    matrix_t **lis = (matrix_t **)malloc(sizeof(matrix_t *) * 2);
    
    m = ((A->nrows - 1) * (A->nrows)) / 2;
    lis[0] = matrix_alloc(m + 1, A->ncols);
    lis[1] = matrix_alloc(m + 1, 1);
    
    for (i = 0; i < A->ncols; ++i)
    {
        lis[0]->data[0][i] = 1;
    }
    lis[1]->data[0][0] = 1;
    n = 1;
    for (i = 0; i < A->nrows; ++i)
    {
        for (j = i + 1; j < A->nrows; ++j)
        {
            for (k = 0; k < A->ncols; ++k)
            {
                lis[0]->data[n][k] = A->data[i][k] - A->data[j][k];
            }
            ++n;
        }
    }
    
    return lis;
}

/* Frees memory used by the Linear system */
void free_lis(matrix_t **lis)
{
    matrix_free(lis[0]);
    matrix_free(lis[1]);
    free(lis);
}

/* Returns 1 if the two support profiles are the same otherwise return 0 */
int is_correct_support(matrix_t *s, matrix_t *supp)
{
    int i;

    for (i = 0; i < s->nrows; ++i)
    {
        //if ((mpq_get_d(supp->data[i][0]) == 0) && (mpq_get_d(s->data[i][0]) != 0))
        if ((supp->data[i][0] == 0) && (s->data[i][0] != 0))
            return 0;
        //if ((mpq_get_d(supp->data[i][0]) == 1) && (mpq_get_d(s->data[i][0]) == 0))
        if ((supp->data[i][0] == 1) && (s->data[i][0] == 0))
            return 0;
    }

    return 1;
}

/* 
 * Given a support 'sol', and the solution to the lis 'val', returns the full
 * vector containing the results of val in their correct positions.
 */
matrix_t *supp_to_sol(matrix_t *sol, matrix_t *val)
{
    int i, k;
    matrix_t *res = matrix_alloc(sol->nrows, 1);

    k = 0;
    for (i = 0; i < sol->nrows; ++i)
        //if (mpq_sgn(sol->data[i][0]) > 0){
        if (sol->data[i][0] > 0) {
            //mpq_set(res->data[i][0], val->data[k][0]);
            res->data[i][0] = val->data[k][0];
            k++;
        }

    return res;
}

/* Check if vector p is indeed a probability vector. */
int check_probability(matrix_t *p)
{
    int res = 0;
    double sum = matrix_sum(p);
    if (err_eq(sum, 1))
        res = 1;

    return res;
}

/* 
 * If current support is a nash equilibrium returns a non-zero value
 * otherwise returns 0.
 */
int is_nash_support()
{
    int err;
    matrix_t *Ap, *Bp, *Bpt, *tmp;
    matrix_t **A, **B;
    
    if (y != NULL) {
        matrix_free(y);
        y = NULL;
    }

    tmp = matrix_row_sub(R, supp_x);
    Ap = matrix_col_sub(tmp, supp_y);
    A = form_lis(Ap);
    matrix_t *yp = lis_solve_system(A[0], A[1], LIS_LU, &err);

    free_lis(A);
    matrix_free(tmp);
    matrix_free(Ap);
    
    if (err || matrix_has_neg(yp))
    {
        err = 1;
    }

    matrix_t *xp;
    if (!err)
    {
        if (x != NULL) {
            matrix_free(x);
            x = NULL;
        }
        Bp = matrix_trans(C);
        tmp = matrix_row_sub(Bp, supp_y);
        Bpt = matrix_col_sub(tmp, supp_x);
        B = form_lis(Bpt);
        xp = lis_solve_system(B[0], B[1], LIS_LU, &err);

        free_lis(B);
        matrix_free(tmp);
        matrix_free(Bp);
        matrix_free(Bpt);
    
        if(err || matrix_has_neg(xp)){
            err = 1;
        }
        if (!err)
            x = supp_to_sol(supp_x, xp);
        matrix_free(xp);
    }

    if(!err){
        y = supp_to_sol(supp_y, yp);
        double eps = compute_epsilon(R, C, x, y, NULL, NULL);
        if (eps < -eps_err || eps > eps_err)
            err = 1;
    }

    matrix_free(yp);

    return !err;
}
