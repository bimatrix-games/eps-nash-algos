#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "lis.h"
#include "checker.h"
#include "../utils/io.h"

#define MAXSTR 100
#define TIMEOUT 15*60

int is_supp;
int is_strategy;
int m;
int n;
int xp, yp;
int is_pure = 1;

/* Allocates memory for the bimatrix */
static void create_vectors()
{
    x = matrix_alloc(m, 1);
    y = matrix_alloc(n, 1);
    supp_x = matrix_alloc(m, 1);
    supp_y = matrix_alloc(n, 1);
}

/* Frees memory used to store the bimatrix */
static void free_vectors()
{
    matrix_free(x);
    matrix_free(y);
    matrix_free(supp_x);
    matrix_free(supp_y);
}

/* Compute the regret when using strategies x, y */
static double check_eps(matrix_t *R, matrix_t *C, int x, int y)
{
    int i, j;
    double tmp;

    tmp = 0;
    for (i = 0; i < m; ++i)
        tmp = fmax(R->data[i][y] - R->data[x][y], tmp);

    for (j = 0; j < n; ++j)
        tmp = fmax(C->data[x][j] - C->data[x][y], tmp);
    return tmp;
}

/* Compute the best pure eps-NE and store in x, y */
static double check_eps_pure(matrix_t *R, matrix_t *C, int *x, int *y)
{
    int i, j, xp, yp;

    xp = 0; yp = 0;
    double eps = 1;
    double tmp;
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

/* Set player p's support to the first n strategies */
static void reset_player(int p, int n)
{
    matrix_t *s = (p == 0) ? supp_x : supp_y;

    int i;
    for (i = 0; i < n; ++i)
        s->data[i][0] = 1;

    for (i = n; i < s->nrows; ++i)
        s->data[i][0] = 0;
}

/* Determine if all of player n's supports have been iterated over. */
static int should_stop(int p, int n)
{
    int i;
    matrix_t *s = (p == 0) ? supp_x : supp_y;

    for (i = s->nrows - n; i < s->nrows; ++i)
        if((s->data[i][0]) == 0)
            return 0;

    return 1;
}

/* Set player p's support to be the next in the binary representation sequence */
static int mat_next(int p)
{
    matrix_t *s = (p == 0) ? supp_x : supp_y;

    int i, carry = 1;
    for (i = 0; i < s->nrows; ++i) {
        int val = (int)(s->data[i][0]);
        val += carry;
        carry = val / 2;
        val %= 2;
        s->data[i][0] = val;
        if (carry == 0)
            break;
    }
    return 1;
}

/*
 * If all of player p's support have been used, then return 0 otherwise
 * set p's support to be the next in the sequence, and return 1
 */
static int player_next(int p, int n)
{
    matrix_t *s = (p == 0) ? supp_x : supp_y;
    if (should_stop(p, n))
        return 0;
    mat_next(p);

    while(matrix_sum(s) != n)
        mat_next(p);

    return 1;
}

/* 
 * Checks all support of size n and returns a 1 if a Nash equilibrium
 * was found otherwise return 0.
 */
static int check_supp(int n)
{
    long int count;
    
    count = 0;
    reset_player(0, n);
    do{
        reset_player(1, n);
        do{
            if (is_nash_support())
                return 1;
        }while (player_next(1, n));
    }while (player_next(0, n));
    printf("\n");
    
    return 0;
}

int supp_size = 1;
/*
 * Computes a Nash equilbrium using support enumeration and returns 1 if a Nash
 * is found otherwise return 0. Due to floating point errors, it is possible that
 * no equilbrium is found on a given game. 
 */
static int supp_enum()
{
    double eps = check_eps_pure(R, C, &xp, &yp);
    if (err_eq(eps, 0)){
        is_pure = 0;
        return 1;
    }
    supp_size = 2;
    int max_size;
    max_size = fmin(R->nrows, R->ncols);
    
    while((supp_size <= max_size) && !check_supp(supp_size))
        ++supp_size;

    if (supp_size > max_size)
        return 0;

    return 1;
}

/* Outputs the current solution to screen */
static void sol_print(matrix_t *x, matrix_t *y)
{    
    int j;
    int count = 0;

    for (j = 0; j < x->nrows; j++)
        if((x->data[j][0]) > 0)
            count++;
    for (j = 0; j < y->nrows; j++)
        if((y->data[j][0]) > 0)
            count++;

    printf("s= %d ", count);
    for (j = 0; j < x->nrows; j++)
        if((x->data[j][0]) > 0)
            printf("%d %.20le ", j+1, (x->data[j][0]));

    for (j = 0; j < y->nrows; j++)
        if((y->data[j][0]) > 0)
            printf("%d %.20le ", j+x->nrows+1, (y->data[j][0]));

    printf("\n");
}

/* 
 * ./supp_enum [i:] 
 * Flags:
 * - '-i file' : Reads the bimatrix game from a file 'filename' using the 
 *   gambit file format otherwise reads in bimatrix from stdin.
 */
int main(int argc, char **argv)
{
    int err = 0;
    int c;

    is_strategy = 1;
    is_supp = 0;
    char *nfg_file = NULL;
    matrix_t **bimatrix;

    while ((c = getopt(argc, argv, "i:ws")) != -1)
        switch (c)
        {
            case 'i':
                nfg_file = optarg;
                break;
            case '?':
                if (isprint(optopt))
                    fprintf(stderr, "Unknown option '-%c.\n", optopt);
                else
                    fprintf(stderr, "Unknown option character '\\x%x.\n",
                            optopt);
                return 1;
            default:
                break;
        }

    if(nfg_file == NULL) {
        bimatrix = read_bimatrix();
        m = bimatrix[0]->nrows;
        n = bimatrix[0]->ncols;
    }
    else {
        FILE *f = fopen(nfg_file, "r");
        bimatrix = read_bimatrix_from_file(f, &m, &n);
        fclose(f);
    }

    create_vectors();
    R = matrix_norm(bimatrix[0]);
    C = matrix_norm(bimatrix[1]);

    int is_nash = supp_enum();

    if (is_pure == 0)
        printf("s= 2 %d 1 %d 1\n", xp + 1, yp + 1 + R->nrows);
    else if (is_nash)
        sol_print(x, y);
    else
        printf("No Nash Found\n");
    
    free_vectors();
    free_matrices(bimatrix);

    return err;
}
