#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "../utils/matrix.h"

int m, n;
matrix_t *R, *C;

/* Read bimatrix from file using gambit file format */
void read_from_file(FILE *f) 
{
	char *buf = (char *) malloc(100 * sizeof(char));
	char c;
	int i, j, tmpn;
	size_t num_bytes = 100;
	double a, b;

	/*
	This checks if we are parsing the correct file type.
	*/
	fgets(buf,num_bytes,f);
	if( strncmp(buf,"NFG 1 D",(size_t)7) != 0 ) {
		fprintf(stderr,"NFG file corrupted, aborting\n");
		exit(1);
	}

	/*
	First we need to ignore all comments (and player names) from the .nfg file.
	*/

	tmpn = 0;
	while(tmpn < 2) { 
		c = fgetc(f);
		if( c == '\"' )
			tmpn++;
	}
	tmpn = 0;
	while(tmpn < 2) {
		c = fgetc(f);
		if( c == '{' || c == '}')
			tmpn++;
	}

	fscanf(f,"%d %d", &m, &n);

	R = matrix_alloc(m, n);
	C = matrix_alloc(m, n);

	fgetc(f); fgetc(f);

	for(i = 0; i < n; i++) {
		for(j = 0; j < m; j++) {
			fscanf(f,"%lf %lf ", &a, &b);
            R->data[j][i] = a;
            C->data[j][i] = b;
		}
	}
	free(buf);
}

/* Compute the regret when using strategies x, y */
double check_eps(matrix_t *R, matrix_t *C, int x, int y)
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
double check_eps_pure(matrix_t *R, matrix_t *C, int *x, int *y)
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

int main(int argc, char **argv)
{
	char c;
	FILE *input;
	
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
	
	read_from_file(input);

    matrix_t *A = matrix_norm(R);
    matrix_t *B = matrix_norm(C);
    matrix_free(R);
    matrix_free(C);
    R = A;
    C = B;

    int x, y;
    check_eps_pure(R, C, &x, &y);

    printf("s= 2 %d 1 %d 1\n", x + 1, y + 1 + m);

	fclose(input);
    matrix_free(R);
    matrix_free(C);
	
	return 0;
}
