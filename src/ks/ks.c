#include <math.h>
#include <float.h>
#include "ks.h"
#include "../utils/util.h"

/* Form LP representing zero-sum game (D, -D) */
cplp_t *form_cplp(matrix_t *D)
{
    int i;
    matrix_t *p = matrix_alloc(D->nrows, 1);
    for (i = 0; i < D->nrows; ++i)
        p->data[i][0] = -1;

    matrix_t *Dp = matrix_augment(D, p);
    matrix_t *s = matrix_alloc(1, Dp->ncols);
    for (i = 0; i < s->ncols - 1; ++i)
        s->data[0][i] = 1;

    matrix_t *Dps = matrix_augment_row(Dp, s);
    matrix_t *b = matrix_alloc(Dps->nrows, 1);
    b->data[b->nrows-1][0] = 1;
    matrix_t *c = matrix_alloc(Dps->ncols, 1);
    c->data[c->nrows - 1][0] = 1;

    cplp_t *lp_t = matrix_to_lp(Dps, b, c, 'L');
    lp_t->obj_sense = CPX_MIN;
    lp_t->sense[lp_t->nrows - 1] = 'E';
    lp_t->lb[Dp->ncols - 1] = -CPX_INFBOUND;

    matrix_free(p);
    matrix_free(Dp);
    matrix_free(s);
    matrix_free(Dps);
    matrix_free(b);
    matrix_free(c);

    return lp_t;
}

/* Create and solve zero-sum game (D, -D) */
cplp_sol_t *solve(matrix_t *D)
{
    cplp_t *cplp = form_cplp(D);
    cplp_sol_t *sol = cplp_solve(cplp);
    cplp_free(cplp);

    return sol;
}

/* Compute the regret when using pure strategies x, y */
double check_eps(matrix_t *R, matrix_t *C, int x, int y)
{
    int i, j;
    double tmp;
    int m = R->nrows;
    int n = R->ncols;

    tmp = 0;
    for (i = 0; i < m; ++i)
        tmp = fmax(R->data[i][y] - R->data[x][y], tmp);

    for (j = 0; j < n; ++j)
        tmp = fmax(C->data[x][j] - C->data[x][y], tmp);
    return tmp;
}

/* Compute the best pure eps-NE and store in x, y */
double check_eps_pure(matrix_t *R, matrix_t *C, double thresh, int *x, int *y)
{
    int i, j, xp, yp;

    xp = 0; yp = 0;
    double eps = 1;
    double tmp;
    int m = R->nrows;
    int n = R->ncols;
    for (i = 0; i < m; ++i){
        if (eps == 0)
            break;
        for (j = 0; j < n; ++j) {
            if (eps == 0)
                break;
            tmp = check_eps(R, C, i, j);
            if (tmp < eps) {
                eps = tmp;
                xp = i;
                yp = j;
            }
        }
    }

    if (x != NULL)
        *x = xp;

    if (y != NULL)
        *y = yp;

    return eps;
}

/* Prints the solution to screen */
void sol_print(sol_t sol)
{    
    int j;
    int count = 0;

    for (j = 0; j < sol.x->nrows; j++)
        if(sol.x->data[j][0] > 0)
            count++;
    for (j = 0; j < sol.y->nrows; j++)
        if(sol.y->data[j][0] > 0)
            count++;

    printf("s= %d ", count);
    for (j = 0; j < sol.x->nrows; j++)
        if(sol.x->data[j][0] > 0)
            printf("%d %.20le ", j+1, sol.x->data[j][0]);

    for (j = 0; j < sol.y->nrows; j++)
        if(sol.y->data[j][0] > 0)
            printf("%d %.20le ", j+sol.x->nrows+1, sol.y->data[j][0]);
    printf("\n");
}

/* Solve the bimatrix game (R, C) using the KS algorithm */
void run_ks(matrix_t *R, matrix_t *C)
{
    int j, x, y;
    double eps1, eps2, thresh;
    sol_t s;
    thresh = (double) 2 / 3;

    matrix_t *A = matrix_norm(R);
    matrix_t *B = matrix_norm(C);
    matrix_t *Dp, *D;
    cplp_sol_t *sol;

    int is_pure = 1;
    eps1 = eps2 = 1;

    eps1 = check_eps_pure(A, B, thresh, &x, &y);
    if (err_le(eps1, 0)) {
        printf("s= 2 %d 1 %d 1 \n", x+1, y+A->nrows+1);
        matrix_free(A);
        matrix_free(B);
        return;
    }

    Dp = matrix_sub(A, B);
    D = matrix_mul_const(Dp, 0.5);
    sol = solve(D);
    is_pure = 0;

    s.x = matrix_alloc(A->nrows, 1);
    s.y = matrix_alloc(A->ncols, 1);
    for (j = 0; j < sol->nrows - 1; j++)
        if(sol->res[1][j] < 0)
            s.x->data[j][0] = -sol->res[1][j];

    for (j = 0; j < sol->ncols - 1; j++)
        if(sol->res[0][j] > 0)
            s.y->data[j][0] = sol->res[0][j];

    cplp_sol_free(sol);

    eps2 = compute_epsilon_supp(A, B, s.x, s.y,NULL, NULL);
    printf("%lf %lf\n", eps1, eps2);
    if (eps1 < eps2)
        printf("s= 2 %d 1 %d 1 \n", x+1, y+A->nrows+1);
    else
        sol_print(s);

    matrix_free(s.x);
    matrix_free(s.y);
    matrix_free(Dp);
    matrix_free(D);

    matrix_free(A);
    matrix_free(B);
}
