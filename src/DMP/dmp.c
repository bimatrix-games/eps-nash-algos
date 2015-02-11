#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

/* Formatting players strategies for output */
#define P1(x) x + 1
#define P2(x) x + m + 1
#define P(x, i) (i == 0) ? P1(x) : P2(x)

int m, n;
double *R, *C;

/* Reads in the bimatrix from file in gambit format */
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

	R = (double*)malloc(m * n * sizeof(double));
	C = (double*)malloc(m * n * sizeof(double));
	fgetc(f); fgetc(f);

	for(i = 0; i < n; i++) {
		for(j = 0; j < m; j++) {
			fscanf(f,"%lf %lf ", &a, &b);
			R[n*j + i] = a;
			C[n*j + i] = b;
		}
	}
	free(buf);
}

/* Gets the best response against player p's strategy i */
int best_response(int i, int p)
{
	int j, b = 0;
	
	if (p == 0)
		for (j = 1; j < n; ++j)
			b = (C[i*n + b] > C[i*n + j]) ? b : j;
	
	else
		for (j = 1; j < m; ++j)
			b = (R[b*n + i] > R[j*n + i]) ? b : j;
	
	return b;
}

/*
 * i and k represent the strategies being used by player p while j is the other
 * players strategy
 */
static int i, j, k, p;

/* Runs the DMP algorithm */
void run_dmp()
{
    if (i <= 0 || i > (m + n))
        i = (rand() % (m + n)) + 1;

    i -= 1;

    p = (i < m) ? 0 : 1;
    i %= m;

	j = best_response(i, p);
	k = best_response(j, (p + 1) % 2);
}

/*
 * ./dmp [i:s:k:]
 *
 * Flags:
 * - '-i file': read in bimatrix game from file
 * - '-s seed': sets the random seed to be used
 * - '-k i': sets the initial strategy of DMP to i
 */
int main(int argc, char **argv)
{
	char c;
	FILE *input;
	
    i = -1;
	while((c = getopt(argc, argv,"i:s:k:")) != -1){
		switch(c) {
		case 'i':
			input = fopen(optarg, "r");
			break;
		case 's':
			srand(atoi(optarg));
			break;
        case 'k':
            i = atoi(optarg);
            break;
		}
	}
	
	read_from_file(input);
	run_dmp();

	if (i == k)
		printf("s= 2 %d 1 %d 1", P(i, p), P(i, (p + 1)%2));
	else
		printf("s= 3 %d 0.5 %d 0.5 %d 1", P(i, p), P(k, p), P(j, (p + 1) % 2));
	printf("\n");

	fclose(input);
	
	return 0;
}
