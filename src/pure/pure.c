#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../utils/io.h"

/* Compute the regret when using strategies x, y */
double check_eps(matrix_t *R, matrix_t *C, int x, int y)
{
    int i, j;
    double tmp;

    tmp = 0;
    for (i = 0; i < R->nrows; ++i)
        tmp = fmax(R->data[i][y] - R->data[x][y], tmp);

    for (j = 0; j < R->ncols; ++j)
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
    for (i = 0; i < R->nrows; ++i){
        if (eps == 0)
            break;
        for (j = 0; j < R->ncols; ++j) {
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

int main(int argc, char **argv)
{
	char c;
    int m, n;
	FILE *input;
    matrix_t **bimatrix;
	
	while((c = getopt(argc, argv,"i:s:k:")) != -1){
		switch(c) {
		case 'i':
            input = fopen(optarg, "r");
			break;
		case 's':
			srand(atoi(optarg));
			break;
        case 'k':
            break;
		}
	}
	
	bimatrix = read_bimatrix_from_file(input, &m, &n);

    matrix_t *R = matrix_norm(bimatrix[0]);
    matrix_t *C = matrix_norm(bimatrix[1]);

    int x, y;
    check_eps_pure(R, C, &x, &y);

    printf("s= 2 %d 1 %d 1\n", x + 1, y + 1 + m);

	fclose(input);
    free_matrices(bimatrix);
    matrix_free(R);
    matrix_free(C);
	
	return 0;
}
