#include "lis.h"

/* Performs forward substitution on Ly = Pb and returns y */
matrix_t* lis_forward(matrix_t *L, matrix_t *Pb)
{
    matrix_t *y = matrix_alloc(Pb->nrows, Pb->ncols);
    double sum, tmp;
    
    int i, j, k;
    
    for (j = 0; j < y->ncols; ++j)
    {
        if(!err_eq(L->data[0][j], 0))
            y->data[0][j] = Pb->data[0][j] / L->data[0][0];
        for (i = 1; i < y->nrows; ++i)
        {
            sum = 0;
            for (k = 0; k < i; ++k)
            {
                tmp = L->data[i][k] * y->data[k][j];
                sum += tmp;
            }
            tmp = Pb->data[i][j] - sum;
            if((L->data[i][i]))
                y->data[i][j] = tmp / L->data[i][i];
        }
    }

    return y;
}

/* Performs backward substitution on Ux = y and returns x */
matrix_t* lis_backward(matrix_t *U, matrix_t *y)
{
    matrix_t *x = matrix_alloc(y->nrows, y->ncols);
    double sum, tmp;

    int i, j, k, l, m, n;
    m = U->nrows;
    n = U->ncols;
    
    for (j = 0; j < x->ncols; ++j) {
        if(!err_eq(U->data[m - 1][m - 1], 0))
            x->data[m - 1][j] = y->data[m - 1][j] / U->data[m - 1][m - 1];

        for (i = m - 2; i >= 0; --i) {
            sum = 0;

            l = i;
            while ((U->data[i][l]) == 0 && l < U->ncols)
                ++l;

            for (k = l + 1; k < m; ++k) {
                tmp = U->data[i][k] * x->data[k][j];
                sum += tmp;
            }

            tmp = y->data[i][j] - sum;
            if(!err_eq(U->data[i][l], 0))
                x->data[l][j] = tmp / U->data[i][l];
        }
    }
    
    return x;
}

/* Solve linear system Ax = b using LU decomposition and returns x */
matrix_t* lis_solve_LU(matrix_t *A, matrix_t *b)
{
    matrix_t *L, *p, *P, *U, *x, *y, *Pb;
    
    L = matrix_alloc(A->nrows, A->nrows);
    p = matrix_alloc(A->nrows, 1);
    U = matrix_alloc(A->nrows, A->ncols);
    
    matrix_LU(A, p, L, U);
    P = matrix_permute(p);
    Pb = matrix_mul(P, b);
    y = lis_forward(L, Pb);
    x = lis_backward(U, y);

    matrix_free(L);
    matrix_free(p);
    matrix_free(U);
    matrix_free(Pb);
    matrix_free(P);
    matrix_free(y);
    
    return x;
}


/* Solves an overdetermined linear system Ax = b using LU decompostition */
matrix_t* lis_solve_over_LU(matrix_t *A, matrix_t *b)
{
    matrix_t *At, *AtA, *Atb, *x;
    
    At = matrix_trans(A);
    AtA = matrix_mul(At, A);
    Atb = matrix_mul(At, b);
    
    x = lis_solve_LU(AtA, Atb);
    
    matrix_free(At);
    matrix_free(AtA);
    matrix_free(Atb);

    return x;
}

/* Solves an underdetermined linear system Ax = b using LU decompostition */
matrix_t* lis_solve_under_LU(matrix_t *A, matrix_t *b)
{
    matrix_t *Ap, *bp, *x;
    
    int i, j, m;
    m = A->ncols;
    Ap = matrix_alloc(m, m);
    bp = matrix_alloc(m, 1);
    
    for (i = 0; i < A->nrows; ++i) {
        bp->data[i][0] = b->data[i][0];
        for (j = 0; j < A->ncols; ++j) 
            Ap->data[i][j] = A->data[i][j];
    }
    
    x = lis_solve_LU(Ap, bp);
    
    matrix_free(Ap);
    matrix_free(bp);

    return x;
}

/* Solves a linear system Ax = b using LU decompostition */
matrix_t* lis_solve_system_LU(matrix_t *A, matrix_t *b, int *err)
{
    matrix_t *x, *bp;
    
    if (A->nrows == A->ncols)
        x = lis_solve_LU(A, b);
    else if(A->nrows < A->ncols)
        x = lis_solve_under_LU(A, b);
    else
        x = lis_solve_over_LU(A, b);
    
    
    bp = matrix_mul(A, x);
    *err = !(matrix_is_equal(b, bp));
    
    matrix_free(bp);

    return x;
}

/* 
 * Solves a linear system of equations Ax = b using the specified method,
 * *err would hold the value 1 if an error occurred during computation
 * otherwise 0.
 */
matrix_t* lis_solve_system(matrix_t *A, matrix_t *b, int method, int *err)
{
    matrix_t *res;
    
    switch (method) {
        case LIS_RREF:
            res = NULL;
            break;
        case LIS_LU:
            res = lis_solve_system_LU(A, b, err);
            break;
    }
    
    return res;
}
