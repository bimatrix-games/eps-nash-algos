#include "bbm.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

/* Form KS LP */
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

/* Create and solve the KS LP */
bbm_sol_t bbm_solve_lp(matrix_t *D)
{
    bbm_sol_t sol;
    cplp_t *cplp = form_cplp(D);
    cplp_sol_t *cplp_sol = cplp_solve(cplp);
    cplp_free(cplp);

    sol.x = matrix_alloc(D->nrows, 1);
    sol.y = matrix_alloc(D->ncols, 1);

    int i;
    for (i = 0; i < D->nrows; ++i)
        sol.x->data[i][0] = -cplp_sol->res[1][i];
    for (i = 0; i < D->ncols; ++i)
        sol.y->data[i][0] = cplp_sol->res[0][i];

    cplp_sol_free(cplp_sol);

    return sol;
}

/* Solve bimatrix game (R,C) using BBM1 */
bbm_sol_t bbm_solve_1(matrix_t *R, matrix_t *C, matrix_t *D)
{
    bbm_sol_t lp_sol = bbm_solve_lp(D);
    matrix_t *x = lp_sol.x;
    matrix_t *y = lp_sol.y;
    double g1, g2;

    bbm_extra_t *ext = bbm_extra_alloc(R, C, x, y);
    double eps = compute_epsilon(ext, &g1, &g2);

    bbm_sol_t sol;

    if (eps <= (3 - sqrt(5)) / 2) {
        sol.x = x;
        sol.y = y;
        bbm_extra_free(ext);
        return sol;
    }

    matrix_t *r1 = matrix_alloc(x->nrows, x->ncols);
    matrix_t *b2 = matrix_alloc(y->nrows, y->ncols);

    double delta = (1 - eps) / (2 - eps);

    if (g1 > g2) {
        r1->data[ext->br_R_y_index][0] = 1;
        sol.x = r1;
        int b;
        matrix_t *r1t_C = matrix_mul_mat_vec(ext->Ct, r1);
        best_response(r1t_C, &b);
        b2->data[b][0] = 1;
        matrix_t *delta_y = matrix_mul_const(y, (1 - delta));
        matrix_t *delta_b2 = matrix_mul_const(b2, delta);
        sol.y = matrix_add(delta_y, delta_b2);

        matrix_free(lp_sol.x);
        matrix_free(lp_sol.y);
        matrix_free(delta_y);
        matrix_free(delta_b2);
        matrix_free(b2);
    } else {
        b2->data[ext->br_Ct_x_index][0] = 1;
        int r;
        matrix_t *R_b2 = matrix_mul_mat_vec(R, b2);
        best_response(R_b2, &r);
        r1->data[r][0] = 1;
        sol.y = b2;
        matrix_t *delta_x = matrix_mul_const(x, (1 - delta));
        matrix_t *delta_r1 = matrix_mul_const(r1, delta);
        sol.x = matrix_add(delta_x, delta_r1);

        matrix_free(lp_sol.x);
        matrix_free(lp_sol.y);
        matrix_free(delta_x);
        matrix_free(delta_r1);
        matrix_free(r1);
    }

    bbm_extra_free(ext);
    matrix_free(x);
    matrix_free(y);
    return sol;
}

/* Compute delta_1 using regret g */
double bbm_delta1(double g)
{
    if (g <= 1.0/3.0 )
        return 0;

    if (g > 0.445043)
        return 1;

    double t1 = 1 / (1 - (2 * g));
    return (1 - g) * (-1 + sqrt(1 + t1 - (1 / g)));
}

/* Compute delta_2 using regret g, delta_1, and h */
double bbm_delta2(double g, double d1, double h)
{
    if (g <= 1.0/3.0)
        return 0;
    if (g > 0.445043)
        return (1 - g) / (2 - g);

    return fmax(0.0, (d1 - g + ((1 - d1) * h)) / (1 + d1 - g));
}

/* 
 * Perform extra step of BBM2 for player p for bimatrix game (R, C) where x is player p's stratey
 * and y is the opponent's strategy.
 */
bbm_sol_t bbm_s(matrix_t *R, matrix_t *C, bbm_extra_t *ext, matrix_t *x, matrix_t *y, double g, int p)
{
    bbm_sol_t sol;

    matrix_t *r1 = matrix_alloc(x->nrows, x->ncols);
    matrix_t *b2 = matrix_alloc(y->nrows, y->ncols);

    if (p == 0) {
        r1->data[ext->br_R_y_index][0] = 1;
    } else {
        r1->data[ext->br_Ct_x_index][0] = 1;
    }

    double delta1 = bbm_delta1(g);

    matrix_t *delta_x = matrix_mul_const(x, (1 - delta1));
    matrix_t *delta_r1 = matrix_mul_const(r1, delta1);
    sol.x = matrix_add(delta_x, delta_r1);

    matrix_t *tmp;
    double h2;
    if (p == 0){
        int b;
        matrix_t *r1t_C = matrix_mul_mat_vec(ext->Ct, sol.x);
        best_response(r1t_C, &b);
        b2->data[b][0] = 1;
        tmp = matrix_mul_mat_vec(ext->xt_C, b2);
        h2 = tmp->data[0][0] - ext->xt_C_y->data[0][0];
    } else {
        int b;
        matrix_t *R_b2 = matrix_mul_mat_vec(R, sol.x);
        best_response(R_b2, &b);
        b2->data[b][0] = 1;
        matrix_t *b2t = matrix_trans(b2);
        tmp = matrix_mul_vec_mat(b2t, ext->R_y);
        h2 = tmp->data[0][0] - ext->xt_R_y->data[0][0];
        matrix_free(b2t);
    }
    double delta2 = bbm_delta2(g, delta1, h2);

    matrix_t *delta_y = matrix_mul_const(y, (1 - delta2));
    matrix_t *delta_b2 = matrix_mul_const(b2, delta2);
    sol.y = matrix_add(delta_y, delta_b2);

    matrix_free(r1);
    matrix_free(b2);
    matrix_free(tmp);
    matrix_free(delta_x);
    matrix_free(delta_r1);
    matrix_free(delta_y);
    matrix_free(delta_b2);

    return sol;
}

/* Solve bimatrix game (R,C) using BBM2 */
bbm_sol_t bbm_solve_2(matrix_t *R, matrix_t *C, matrix_t *D)
{
    bbm_sol_t lp_sol = bbm_solve_lp(D);
    matrix_t *x = lp_sol.x;
    matrix_t *y = lp_sol.y;

    double g1, g2;

    bbm_extra_t *ext = bbm_extra_alloc(R, C, x, y);
    double eps = compute_epsilon(ext, &g1, &g2);
    double thresh = (double)1/3;

    if (eps <= thresh)
        return lp_sol;

    bbm_sol_t sol;
    if (g1 >= g2) {
        sol = bbm_s(R, C, ext, x, y, g1, 0);
    }
    else {
        sol = bbm_s(R, C, ext, y, x, g2, 1);
        matrix_t *t = sol.x;
        sol.x = sol.y;
        sol.y = t;
    }

    matrix_free(x);
    matrix_free(y);
    return sol;
}

/* Output the solution to screen */
void bbm_sol_print(bbm_sol_t sol)
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

/* Run bbm */
void run_bbm(matrix_t *R, matrix_t *C, int bbm)
{
    matrix_t *A = matrix_norm(R);
    matrix_t *B = matrix_norm(C);

    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    matrix_t *Dp = matrix_sub(A, B);
    matrix_t *D = matrix_mul_const(Dp, 0.5);

    bbm_sol_t sol;
    if (bbm == 0)
        sol = bbm_solve_1(A, B, D);
    else
        sol = bbm_solve_2(A, B, D);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;

    printf("Compute Time %lf\n", elapsed);
    bbm_sol_print(sol);

    matrix_free(sol.x);
    matrix_free(sol.y);
    matrix_free(A);
    matrix_free(B);
    matrix_free(Dp);
    matrix_free(D);
}
