#include "ks.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "../utils/util.h"

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

/* Form LP representing zero-sum game (D, -D) */
cplp_t *form_ks_cplp(matrix_t *D)
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
    matrix_free(Dps);
    matrix_free(s);
    matrix_free(c);
    matrix_free(b);

    return lp_t;
}

/* Create and solve zero-sum game (D, -D) */
sol_t solve_ks(matrix_t *D)
{
    sol_t sol;
    cplp_t *cplp = form_ks_cplp(D);
    cplp_sol_t *lp_sol = cplp_solve(cplp);
    cplp_free(cplp);

    sol.x = matrix_alloc(D->nrows, 1);
    sol.y = matrix_alloc(D->ncols, 1);
    int i;
    for (i = 0; i < D->nrows; ++i)
        sol.x->data[i][0] = -lp_sol->res[1][i];
    for (i = 0; i < D->ncols; ++i)
        sol.y->data[i][0] = lp_sol->res[0][i];
    cplp_sol_free(lp_sol);

    return sol;
}

/* 
 * Forms the LP to find the best WSNE with a given pair of supports Sr and Sc
 * for the row and column players respectively.
 */
cplp_t *form_imp_lp(matrix_t *R, matrix_t *Sr, matrix_t *Sc)
{
    int i, ip, j, k, m, l;
    m = (R->nrows * Sr->nrows) - Sr->nrows;
    matrix_t *M = matrix_alloc(m+2, R->ncols + 1);
    matrix_t *b = matrix_alloc(m+2, 1);
    matrix_t *c = matrix_alloc(R->ncols+1, 1);
    c->data[c->nrows - 1][0] = 1;

    b->data[m][0] = 1;
    b->data[m+1][0] = 1;

    k = 0;
    for (i = 0; i < Sr->nrows; ++i) {
        ip = Sr->data[i][0];
        for (j = 0; j < R->nrows; ++j) {
            if (j == ip)
                continue;
            for (l = 0; l <  R->ncols; ++l){
                M->data[k][l] = R->data[j][l] - R->data[ip][l];
            }
            M->data[k][l] = -1;
            k++;
        }
    }

    k = 0;
    for (i = 0; i < R->ncols; ++i){
        M->data[m][i] = 1;
        if (k < Sc->nrows && i == Sc->data[k][0]) {
            M->data[m+1][i] = 1;
            k++;
        }
    }

    cplp_t *lp_t = matrix_to_lp(M, b, c, 'L');
    lp_t->obj_sense = CPX_MIN;
    lp_t->sense[lp_t->nrows - 1] = 'E';
    lp_t->sense[lp_t->nrows - 2] = 'E';

    matrix_free(M);
    matrix_free(b);
    matrix_free(c);
    return lp_t;
}

/* 
 * Given a set of supports Sr and Sc for the row and column players respectively,
 * returns the best WSNE possible with those set of supports.
 */
sol_t solve_imp(matrix_t *R, matrix_t *C, matrix_t *Ct, matrix_t *Sr, matrix_t *Sc, double *eps)
{
    sol_t sol;
    sol.x = matrix_alloc(R->nrows, 1);
    sol.y = matrix_alloc(R->ncols, 1);

    cplp_t *imp_lp = form_imp_lp(R, Sr, Sc);
    cplp_sol_t *imp_sol = cplp_solve(imp_lp);
    cplp_free(imp_lp);

    int i;
    for (i = 0; i < R->ncols; ++i)
        sol.y->data[i][0] = imp_sol->res[0][i];
    
    *eps = imp_sol->obj_val;
    cplp_sol_free(imp_sol);

    imp_lp = form_imp_lp(Ct, Sc, Sr);
    imp_sol = cplp_solve(imp_lp);
    for (i = 0; i < R->nrows; ++i)
        sol.x->data[i][0] = imp_sol->res[0][i];
    *eps = fmax(*eps, imp_sol->obj_val);

    cplp_sol_free(imp_sol);
    cplp_free(imp_lp);

    return sol;
}

/* Find the best WSNE with a 2 x 2 support pair */
sol_t solve_two_supp(matrix_t *R, matrix_t *C, matrix_t *Ct, double *e)
{
    sol_t sol = {NULL, NULL};

    double b_eps = 1;
    double eps = 1;
    int i, j, k, l;
    matrix_t *Sr = matrix_alloc(2, 1);
    matrix_t *Sc = matrix_alloc(2, 1);

    for (i = 0; i < R->nrows - 1; ++i) {
        for (j = i + 1; j < R->nrows; ++j) {
            for (k = 0; k < R->ncols - 1; ++k) {
                for (l = k + 1; l < R->ncols; ++l) {
                    Sr->data[0][0] = i;
                    Sr->data[1][0] = j;
                    Sc->data[0][0] = k;
                    Sc->data[1][0] = l;
                    sol_t s = solve_imp(R, C, Ct, Sr, Sc, &eps);
                    if (eps < b_eps){
                        if (sol.x != NULL) {
                            matrix_free(sol.x);
                            matrix_free(sol.y);
                        }
                        sol = s;
                        b_eps = eps;
                    }
                    else{
                        matrix_free(s.x);
                        matrix_free(s.y);
                    }
                }
            }
        }
    }
    *e = b_eps;
    matrix_free(Sr);matrix_free(Sc);

    return sol;
}

/* 
 * Solve the zero-sum game (D, -D) and return the best WSNE with the support
 * from the zero-sum game solution
 */
sol_t solve_ks_imp(matrix_t *R, matrix_t *C, matrix_t *Ct, double *eps)
{
    matrix_t *Dp = matrix_sub(R, C);
    matrix_t *D = matrix_mul_const(Dp, 0.5);
    sol_t sol = solve_ks(D);

    int i, j, k;
    k = 0;
    
    for (i = 0; i < R->nrows; ++i)
        if (sol.x->data[i][0] > 0)
            ++k;

    matrix_t *Sr = matrix_alloc(k, 1);
    j = 0;
    for (i = 0; i < R->nrows && j < k; ++i)
        if (sol.x->data[i][0] > 0)
            Sr->data[j++][0] = i;

    k = 0;
    for (i = 0; i < R->ncols; ++i)
        if (sol.y->data[i][0] > 0)
            ++k;

    matrix_t *Sc = matrix_alloc(k, 1);
    j = 0;
    for (i = 0; i < R->ncols && j < k; ++i)
        if (sol.y->data[i][0] > 0)
            Sc->data[j++][0] = i;

    sol_t s = solve_imp(R, C, Ct, Sr, Sc, eps);

    matrix_free(Sr);
    matrix_free(Sc);
    matrix_free(Dp);
    matrix_free(D);
    matrix_free(sol.x);
    matrix_free(sol.y);
    return s;
}

/* Compute the regret when using pure strategies x, y */
double check_eps(matrix_t *R, matrix_t *C, int x, int y)
{
    int i, j;
    double tmp;

    tmp = 0;
    for (i = 0; i < R->nrows; ++i)
        tmp = fmax(R->data[i][y] - R->data[x][y], tmp);

    for (j = 0; j < C->ncols; ++j)
        tmp = fmax(C->data[x][j] - C->data[x][y], tmp);

    return tmp;
}

/* Compute the best pure eps-NE and store in x, y */
double check_eps_pure(matrix_t *R, matrix_t *C, int *x, int *y)
{
    int i, j, xp, yp;

    xp = 0; yp = 0;
    double eps = 1;
    double tmp;
    for (i = 0; i < R->nrows; ++i)
        for (j = 0; j < R->ncols; ++j) {
            tmp = check_eps(R, C, i, j);
            if (tmp < eps) {
                eps = tmp;
                xp = i;
                yp = j;
            }
        }

    if (x != NULL)
        *x = xp;

    if (y != NULL)
        *y = yp;

    return eps;
}

/* Solve the bimatrix game (R, C) using the KS+ algorithm */
sol_t solve(matrix_t *R, matrix_t *C)
{
    int x, y;
    double eps_pure = check_eps_pure(R, C, &x, &y);
    sol_t sol_pure;
    sol_pure.x = matrix_alloc(R->nrows, 1);
    sol_pure.y = matrix_alloc(R->ncols, 1);
    sol_pure.x->data[x][0] = 1;
    sol_pure.y->data[y][0] = 1;

    matrix_t *Ct = matrix_trans(C);
    double eps_ts, eps_ks_i;
    sol_t sol_two = solve_two_supp(R, C, Ct, &eps_ts);
    sol_t sol_ks_i = solve_ks_imp(R, C, Ct, &eps_ks_i);
    matrix_free(Ct);

    if (eps_pure <= eps_ts && eps_pure <= eps_ks_i){
        matrix_free(sol_two.x);
        matrix_free(sol_two.y);
        matrix_free(sol_ks_i.x);
        matrix_free(sol_ks_i.y);
        return sol_pure;
    }
    else if (eps_ts < eps_pure && eps_ts <= eps_ks_i){
        matrix_free(sol_pure.x);
        matrix_free(sol_pure.y);
        matrix_free(sol_ks_i.x);
        matrix_free(sol_ks_i.y);
        return sol_two;
    }
    else {
        matrix_free(sol_pure.x);
        matrix_free(sol_pure.y);
        matrix_free(sol_two.x);
        matrix_free(sol_two.y);
        return sol_ks_i;
    }
}

/* Solve the bimatrix game (R, C) using the KS+ algorithm */
void run_ks(matrix_t *R, matrix_t *C)
{
    matrix_t *A = matrix_norm(R);
    matrix_t *B = matrix_norm(C);

    sol_t s = solve(A,B);
    sol_print(s);

    matrix_free(A);
    matrix_free(B);
    matrix_free(s.x);
    matrix_free(s.y);
}
