#include "ts.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

/* Prints out the solution to screen */
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

/* Forms LP to equalize players regret */
cplp_t *form_equal_regret_lp(int m, int n, matrix_t *m1, matrix_t *m2, matrix_t *m3)
{
    int i;

    matrix_t *zeros = matrix_alloc(n, 1);
    matrix_t *m_ones = matrix_alloc(n, 1);
    for (i = 0; i < m_ones->nrows; ++i)
        m_ones->data[i][0] = -1;

    matrix_t *zero = matrix_alloc(1, 1);
    matrix_t *one = matrix_alloc(1, 1);
    one->data[0][0] = 1;
    matrix_t *m_one = matrix_alloc(1, 1);
    m_one->data[0][0] = -1;

    matrix_t *m4 = matrix_trans(m1);

    double max;
    matrix_t *max_m1 = matrix_alloc(1, 1);
    matrix_min_max(m1, NULL, &max);
    max_m1->data[0][0] = max;

    matrix_t *C1 = matrix_alloc(1, m + 2);
    C1->data[0][m] = 1;

    matrix_t *C2 = matrix_augment_cols(3, m2, zeros, m_ones);

    matrix_t *C3 = matrix_augment_cols(3, m3, m_one, one);

    matrix_t *C4 = matrix_alloc(1, m + 2);
    for (i = 0; i < m; ++i)
        C4->data[0][i] = 1;

    matrix_t *M = matrix_augment_rows(4, C1, C2, C3, C4);
    matrix_t *rhs = matrix_augment_rows(4, max_m1, zeros, zero, one);
    matrix_t *m_m1 = matrix_mul_const(m1, -1);
    matrix_t *c = matrix_augment_rows(3, m_m1, one, zero);

    cplp_t *lp = matrix_to_lp(M, rhs, c, 'L');
    lp->obj_sense=CPX_MIN;
    lp->sense[0] = 'E';
    lp->sense[lp->nrows - 1] = 'E';

    matrix_free(zeros);
    matrix_free(zero);
    matrix_free(one);
    matrix_free(m_one);
    matrix_free(m_ones);
    matrix_free(m4);
    matrix_free(max_m1);
    matrix_free(C1);
    matrix_free(C2);
    matrix_free(C3);
    matrix_free(C4);
    matrix_free(M);
    matrix_free(rhs);
    matrix_free(m_m1);
    matrix_free(c);
    return lp;
}

/* Equalize the players regret */
sol_t equal_regret(matrix_t *A, matrix_t *B, matrix_t *Bt, matrix_t *x, matrix_t *y, extra_t *ext, int *eq)
{
    *eq = 1;
    double er, ec;
    compute_epsilon_extra(ext, &er, &ec);
    cplp_t *lp;
    sol_t sol;

    if (er == ec) {
        *eq = 0;
        return sol;
    }
    matrix_t *m1, *m2, *m3;
    if (er > ec) {
        m1 = ext->R_y;
        m2 = Bt;
        m3 = matrix_sub(ext->R_y, ext->C_y);
        matrix_t *m3t = matrix_trans(m3);
        lp = form_equal_regret_lp(x->nrows, y->nrows, m1, m2, m3t);
        matrix_free(m3t);
        matrix_free(m3);
    } else {
        m1 = ext->Ct_x;
        m2 = A;
        m3 = matrix_sub(ext->xt_C, ext->xt_R);
        lp = form_equal_regret_lp(y->nrows, x->nrows, m1, m2, m3);
        matrix_free(m3);
    }

    cplp_sol_t *lp_sol = cplp_solve(lp);
    cplp_free(lp);

    int i;
    if (er > ec) {
        sol.x = matrix_alloc(x->nrows, 1);
        for (i = 0; i < x->nrows; ++i)
            sol.x->data[i][0] = lp_sol->res[0][i];
        sol.y = matrix_copy(y);
    } else if(ec > er) {
        sol.x = matrix_copy(x);
        sol.y = matrix_alloc(y->nrows, 1);
        for (i = 0; i < y->nrows; ++i)
            sol.y->data[i][0] = lp_sol->res[0][i];
    } else {
        *eq = 0;
    }
    cplp_sol_free(lp_sol);

    return sol;
}

/* Forms the steepest descent LP */
cplp_t *form_steepest_lp(matrix_t *A, matrix_t *B, matrix_t *Bt, matrix_t *x, matrix_t *y, extra_t *ext, double delta, int *lp_size)
{
    int i;
    matrix_t *ones, *zeros, *tmp;
    matrix_t *r = find_delta_br(ext->R_y, delta);
    matrix_t *s = find_delta_br(ext->Ct_x, delta);
    matrix_t *A0 = matrix_alloc(s->nrows, A->ncols);
    matrix_t *B0 = matrix_alloc(r->nrows, A->nrows);
    
    ones = matrix_alloc(r->nrows, 1);
    for (i = 0; i < r->nrows; ++i)
        ones->data[i][0] = -1;
    zeros = matrix_alloc(r->nrows, 1);
    matrix_t *A_r = matrix_row_sub_index(A, r);
    matrix_t *C1 = matrix_augment_cols(5, B0, A_r, ones, zeros, zeros);
    matrix_free(A_r);
    matrix_free(ones);
    matrix_free(zeros);

    ones = matrix_alloc(s->nrows, 1);
    for (i = 0; i < s->nrows; ++i)
        ones->data[i][0] = -1;
    zeros = matrix_alloc(s->nrows, 1);
    matrix_t *Bt_s = matrix_row_sub_index(Bt, s);
    matrix_t *C2 = matrix_augment_cols(5, Bt_s, A0, zeros, ones, zeros);
    matrix_free(ones);
    matrix_free(zeros);
    matrix_free(Bt_s);

    zeros = matrix_alloc(1, 1);
    ones = matrix_alloc(1, 1);
    ones->data[0][0] = 1;
    matrix_t *n_ones = matrix_alloc(1, 1);
    n_ones->data[0][0] = -1;
    tmp = matrix_trans(ext->R_y);
    matrix_t *t1 = matrix_neg(tmp);
    matrix_t *t2 = matrix_neg(ext->xt_R);
    matrix_t *C3 = matrix_augment_cols(5, t1, t2, ones, zeros, n_ones);
    matrix_free(tmp);
    matrix_free(t1);
    matrix_free(t2);

    tmp = matrix_trans(ext->C_y);
    t1 = matrix_neg(tmp);
    t2 = matrix_neg(ext->xt_C);
    matrix_t *C4 = matrix_augment_cols(5, t1, t2, zeros, ones, n_ones);
    matrix_free(tmp);
    matrix_free(t1);
    matrix_free(t2);
    matrix_free(zeros);
    matrix_free(ones);
    matrix_free(n_ones);

    matrix_t *C5 = matrix_alloc(1, C4->ncols);
    for (i = 0; i < x->nrows; ++i)
        C5->data[0][i] = 1;

    matrix_t *C6 = matrix_alloc(1, C4->ncols);
    for (i = x->nrows; i < x->nrows + y->nrows; ++i)
        C6->data[0][i] = 1;

    matrix_t *Ap = matrix_augment_rows(6, C1, C2, C3, C4, C5, C6);
    matrix_t *b = matrix_alloc(Ap->nrows, 1);
    b->data[Ap->nrows-4][0] = -(ext->xt_R_y->data[0][0]);
    b->data[Ap->nrows-3][0] = -(ext->xt_C_y->data[0][0]);
    b->data[Ap->nrows-2][0] = 1;
    b->data[Ap->nrows-1][0] = 1;
    matrix_t *c = matrix_alloc(Ap->ncols, 1);
    c->data[c->nrows - 1][0] = 1;
    *lp_size = Ap->nrows;

    cplp_t *lp_t = matrix_to_lp(Ap, b, c, 'L');
    lp_t->obj_sense = CPX_MIN;
    lp_t->sense[lp_t->nrows - 1] = 'E';
    lp_t->sense[lp_t->nrows - 2] = 'E';
    lp_t->lb[Ap->ncols - 1] = -CPX_INFBOUND;

    matrix_free(Ap);
    matrix_free(b);
    matrix_free(c);
    matrix_free(A0);
    matrix_free(B0);
    matrix_free(r);
    matrix_free(s);
    matrix_free(C1);
    matrix_free(C2);
    matrix_free(C3);
    matrix_free(C4);
    matrix_free(C5);
    matrix_free(C6);

    return lp_t;
}

/* Computes and returns the steepest descent */
steep_sol_t steepest_descent(matrix_t *A, matrix_t *B, matrix_t *Bt, matrix_t *x, matrix_t *y, extra_t *ext, double delta, int *lp_size)
{
    cplp_t *lp = form_steepest_lp(A, B, Bt, x, y, ext, delta, lp_size);
    cplp_sol_t *lp_sol = cplp_solve(lp);
    cplp_free(lp);

    sol_t prob;
    prob.x = matrix_alloc(A->nrows, 1);
    prob.y = matrix_alloc(A->ncols, 1);
    
    int i;
    for (i = 0; i < A->nrows; ++i)
        prob.x->data[i][0] = lp_sol->res[0][i];
    for (; i < A->nrows + A->ncols; ++i)
        prob.y->data[i - A->nrows][0] = lp_sol->res[0][i];

    steep_sol_t sol = {prob, lp_sol};

    return sol;
}

/* Computes the new value of x* from x in direction xp at distance eps */
matrix_t *new_point_i(double eps, matrix_t *x, matrix_t *xp)
{
    matrix_t *t1, *t2, *res;

    t1 = matrix_sub(xp, x);
    t2 = matrix_mul_const(t1, eps);
    res = matrix_add(x, t2);

    matrix_free(t1);
    matrix_free(t2);
    return res;
}

/* Calculates the new point from (x, y) in direction (xp, yp) at distance eps */
sol_t new_point(double eps, matrix_t *x, matrix_t *xp, matrix_t *y, matrix_t *yp)
{
    sol_t sol;
    sol.x = new_point_i(eps, x, xp);
    sol.y = new_point_i(eps, y, yp);
    return sol;
}

/* 
 * Computes the value of delta_Df_i where i = (p+1), x, xp (y, yp if p = 1) are player 1's
 * strategy and direction respectively, and A is R when p = 0 otherwise A=C^T 
 */
double delta_Df_i(matrix_t *A, matrix_t *x, matrix_t *xp, matrix_t *y, matrix_t *yp, extra_t *ext, double delta, int p)
{
    double xt_A_y, max_Ayp_s;
    matrix_t *Ay, *xt;
    if (p == 0) {
        xt_A_y = ext->xt_R_y->data[0][0];
        Ay = ext->R_y;
    }
    else {
        xt_A_y = ext->xt_C_y->data[0][0];
        Ay = ext->Ct_x;
    }
    xt = matrix_trans(x);

    matrix_t *Ayp = matrix_mul_mat_vec(A, yp);
    matrix_t *S_Ay = find_delta_br(Ay, delta);
    matrix_t *Ayp_S_Ay = matrix_row_sub_index(Ayp, S_Ay);
    matrix_t *xpt = matrix_trans(xp);
    matrix_t *xpt_Ay = matrix_mul_vec_mat(xpt, Ay);
    matrix_t *xt_Ayp = matrix_mul_vec_mat(xt, Ayp);
    matrix_min_max(Ayp_S_Ay, NULL, &max_Ayp_s);

    double res = max_Ayp_s - xt_Ayp->data[0][0] - xpt_Ay->data[0][0] + xt_A_y;

    matrix_free(xt);
    matrix_free(Ayp);
    matrix_free(S_Ay);
    matrix_free(Ayp_S_Ay);
    matrix_free(xpt);
    matrix_free(xpt_Ay);
    matrix_free(xt_Ayp);

    return res;
}

/* Returns the value of the gradient at (x, y) in the direction (xp, yp) */
double delta_Df(matrix_t *A, matrix_t *B,matrix_t *Bt, matrix_t *x, matrix_t *y, matrix_t *xp, matrix_t *yp, extra_t *ext, double delta)
{
    double df1, df2, eps, er, ec;

    eps = compute_epsilon_extra(ext, &er, &ec);
    df1 = delta_Df_i(A, x, xp, y, yp, ext, delta, 0);
    df2 = delta_Df_i(Bt, y, yp, x, xp, ext, delta, 1);

    if (er > ec)
        return df1 - eps;
    else if (ec > er)
        return df2 - eps;

    return fmax(df1,df2) - eps;
}

/* 
 * Given the solution to the steepest descent LP and the indices of rows and columns
 * in the delta best response, extracts and returns the values of w and z from the
 * dual solutions.
 */
sol_t extract_wz_from_dual(steep_sol_t steep, matrix_t *r, matrix_t *s)
{
    sol_t sol;

    matrix_t *x = matrix_alloc(steep.prob.x->nrows, 1);
    matrix_t *y = matrix_alloc(steep.prob.y->nrows, 1);

    int i, j, k;
    j = 0;
    for (i = 0; i < r->nrows; i++) {
        k = r->data[i][0];
        x->data[k][0] = -steep.lp_sol->res[1][j++];
    }

    for (i = 0; i < s->nrows; i++) {
        k = s->data[i][0];
        y->data[k][0] = -steep.lp_sol->res[1][j++];
    }

    sol.x = matrix_prob_norm(x);
    sol.y = matrix_prob_norm(y);
    matrix_free(x);
    matrix_free(y);

    return sol;
}

/* Form the LP needed to compute the value of lambda (mu for player 2) */
cplp_t *form_prime_lp(matrix_t *R, matrix_t *w, matrix_t *x, matrix_t *S_x)
{
    matrix_t *w_x = matrix_sub(w, x);
    matrix_t *wx_t = matrix_trans(w_x);
    matrix_t *wx_t_R = matrix_mul_vec_mat(wx_t, R);
    matrix_t *c = matrix_trans(wx_t_R);

    matrix_t *A = matrix_alloc(2, wx_t_R->ncols);
    matrix_t *b = matrix_alloc(2, 1);
    b->data[0][0] = 1;
    b->data[1][0] = 1;

    int i, j, k;
    j = k = 0;
    for (i = 0; i < A->ncols; ++i) {
        A->data[0][i] = 1;
        if (k < S_x->nrows && i == S_x->data[k][0]) {
            A->data[1][i] = 1;
            k++;
        }
    }

    cplp_t *lp = matrix_to_lp(A, b, c, 'E');
    lp->obj_sense=CPX_MIN;

    matrix_free(A);
    matrix_free(b);
    matrix_free(c);
    matrix_free(wx_t_R);
    matrix_free(wx_t);
    matrix_free(w_x);

    return lp;
}

/* Computes the value of lambda (mu for player 2) */
matrix_t *find_prime(matrix_t *A, matrix_t *w, matrix_t *x, matrix_t *S_x, double *lambda)
{
    cplp_t *lp = form_prime_lp(A, w, x, S_x);
    cplp_sol_t *lp_sol = cplp_solve(lp);
    cplp_free(lp);

    matrix_t *y = matrix_alloc(A->ncols, 1);

    int i;
    for (i = 0; i < A->ncols; ++i)
        y->data[i][0] = lp_sol->res[0][i];

    *lambda = lp_sol->obj_val;
    cplp_sol_free(lp_sol);
    return y;
}

/* Computes the values of lambda and mu */
sol_t find_primes(matrix_t *A, matrix_t *B, matrix_t *Bt, matrix_t *x, matrix_t *y, matrix_t *w, matrix_t *z, extra_t *ext, double delta, double *lambda, double *mu)
{
    sol_t sol;

    matrix_t *S_x = find_delta_br(ext->Ct_x, delta);
    sol.y = find_prime(A, w, x, S_x, lambda);
    matrix_t *S_y = find_delta_br(ext->R_y, delta);
    sol.x = find_prime(Bt, z, y, S_y, mu);

    matrix_free(S_x);
    matrix_free(S_y);

    return sol;
}

/* 
 * Compute the value of x~ (y~ for the second player) given x*, w*, lambda and mu
 * (z*, y*, mu and lambda for player 2)
 */
matrix_t *extra_point_i(matrix_t *x, matrix_t *w, double lambda, double mu)
{
    double w1 = 1 / (1 + lambda - mu);
    double x1 = (lambda - mu) / (1 + lambda - mu);
    matrix_t *wp = matrix_mul_const(w, w1);
    matrix_t *xp = matrix_mul_const(x, x1);
    matrix_t *res = matrix_add(wp, xp);

    matrix_free(wp);
    matrix_free(xp);

    return res;
}

/* 
 * Extract the extra point (x~, y~) using x, y and w, z given the values of the
 * LP's lambda and mu
 */ 
sol_t extra_point(matrix_t *x, matrix_t *y, matrix_t *w, matrix_t *z, double lambda, double mu)
{
    sol_t sol;

    if (lambda >= mu) {
        sol.x = extra_point_i(x, w, lambda, mu);
        sol.y = matrix_copy(z);
    } else {
        sol.x = matrix_copy(w);
        sol.y = extra_point_i(y, z, mu, lambda);
    }

    return sol;
}

/* 
 * Solves the bimatrix game (A, B) using the starting point (x0, y0), with a 
 * fixed delta and returns the final point of the steepest descent and the 
 * extra points.
 */
ts_sol_t solve_ts(matrix_t *A, matrix_t *B, matrix_t *Bt, matrix_t *x0, matrix_t *y0, double delta, long *size)
{
    matrix_t *xp, *yp;
    sol_t np;
    extra_t *ext;
    long count = 0;
    steep_sol_t steep;
    matrix_t *x, *y;
    x = matrix_copy(x0);
    y = matrix_copy(y0);
    double e = delta / (delta + 1);
    while (1) {
        count++;
        ext = extra_alloc(A, B, Bt, x, y);
        int eq;
        sol_t eq_reg = equal_regret(A, B, Bt, x, y, ext, &eq);
        if (eq != 0) {
            matrix_free(x);
            matrix_free(y);
            extra_free(ext);
            x = eq_reg.x;
            y = eq_reg.y;
            ext = extra_alloc(A, B, Bt, x, y);
        }

        int lp_size;
        steep = steepest_descent(A, B, Bt, x, y, ext, delta, &lp_size);
        xp = steep.prob.x;
        yp = steep.prob.y;
        *size += lp_size;

        double Df = delta_Df(A, B, Bt, x, y, xp, yp, ext, delta);
        if (Df >= -delta)
            break;

        np = new_point(e, x, xp, y, yp);
        extra_free(ext);
        cplp_sol_free(steep.lp_sol);
        matrix_free(xp);
        matrix_free(yp);
        matrix_free(x);
        matrix_free(y);

        x = np.x;
        y = np.y;
    }

    matrix_t *r = find_delta_br(ext->R_y, delta);
    matrix_t *s = find_delta_br(ext->Ct_x, delta);
    sol_t wz = extract_wz_from_dual(steep, r, s);
    matrix_free(r);
    matrix_free(s);

    double lambda, mu;
    sol_t xp_yp = find_primes(A, B, Bt, x, y, wz.x, wz.y, ext, delta, &lambda, &mu);
    sol_t ext_point = extra_point(x, y, wz.x, wz.y, lambda, mu);

    extra_free(ext);
    ts_sol_t ts_sol;
    ts_sol.x_y.x = x;
    ts_sol.x_y.y = y;
    ts_sol.xp_yp.x = xp;
    ts_sol.xp_yp.y = yp;
    double e1 = compute_epsilon(A, B, x, y, NULL, NULL);
    double e2 = compute_epsilon(A, B, ext_point.x, ext_point.y, NULL, NULL);
    if (e1 < e2) {
        ts_sol.xh_yh.x = matrix_copy(x);
        ts_sol.xh_yh.y = matrix_copy(y);
        matrix_free(ext_point.x);
        matrix_free(ext_point.y);
    } else {
        ts_sol.xh_yh = ext_point;
    }
    ts_sol.iterations = count;
    cplp_sol_free(steep.lp_sol);
    matrix_free(wz.x);
    matrix_free(wz.y);
    matrix_free(xp_yp.x);
    matrix_free(xp_yp.y);

    return ts_sol;
}

/* Solves the bimatrix game (R, C) using the TS algorithm with a given delta */
void run_ts(matrix_t *R, matrix_t *C, double delta)
{
    /* Normalize the payoff matrices */
    matrix_t *A = matrix_norm(R);
    matrix_t *B = matrix_norm(C);
    matrix_t *Bt = matrix_trans(B);

    /* Initialize the starting point of TS to be a randomly chosen point */
    matrix_t *x = matrix_alloc(A->nrows, 1);
    matrix_t *y = matrix_alloc(A->ncols, 1);
    matrix_rand(x);
    matrix_rand(y);
    matrix_t *tx = matrix_prob_norm(x);
    matrix_t *ty = matrix_prob_norm(y);
    matrix_free(x);
    matrix_free(y);
    x = tx;
    y = ty;

    long size = 0;

    /* Get the result of the TS algorithm starting from point x, y with a
     * fixed delta, storing the sum of the LP sizes of all iterations */
    ts_sol_t sol = solve_ts(A, B, Bt, x, y, delta, &size);
    matrix_free(x);
    matrix_free(y);

    sol_t s;
    double e1, e2, e3, e4, e5;
    /* Computes the best epsilon from the extra points */
    e1 = compute_epsilon(A, B, sol.x_y.x, sol.x_y.y, NULL, NULL);
    e2 = compute_epsilon(A, B, sol.xp_yp.x, sol.xp_yp.y, NULL, NULL);
    e3 = compute_epsilon(A, B, sol.xh_yh.x, sol.xh_yh.y, NULL, NULL);
    e4 = compute_epsilon(A, B, sol.xp_yp.x, sol.x_y.y, NULL, NULL);
    e5 = compute_epsilon(A, B, sol.x_y.x, sol.xp_yp.y, NULL, NULL);

    char *out;
    if (e1 <= e2 && e1 <= e3 && e1 <= e4 && e1 <= e5){
        s = sol.x_y;
        out = "x_y";
    }
    if (e2 < e1 && e2 <= e3 && e2 <= e4 && e2 <= e5){
        s = sol.xp_yp;
        out = "xp_yp";
    }
    if (e3 < e2 && e3 < e1 && e3 <= e4 && e3 <= e5){
        s = sol.xh_yh;
        out = "xh_yh";
    }
    if (e4 < e2 && e4 < e3 && e4 < e1 && e4 <= e5){
        s = (sol_t){.x = sol.xp_yp.x, .y = sol.x_y.y};
        out = "xp_y";
    }
    if (e5 < e2 && e5 < e3 && e5 < e4 && e5 < e1){
        s = (sol_t){.x = sol.x_y.x, .y = sol.xp_yp.y};
        out = "x_yp";
    }

    sol_print(s);
    printf("Iterations %ld\n", sol.iterations);
    printf("Solution %s\n", out);
    printf("Eps %lf %lf %lf %lf %lf\n", e1, e2, e3, e4, e5);
    printf("Avg LP Size %lf\n", (double) size / sol.iterations);
    matrix_free(sol.x_y.x);
    matrix_free(sol.x_y.y);
    matrix_free(sol.xp_yp.x);
    matrix_free(sol.xp_yp.y);
    matrix_free(sol.xh_yh.x);
    matrix_free(sol.xh_yh.y);

    matrix_free(A);
    matrix_free(B);
    matrix_free(Bt);
}
