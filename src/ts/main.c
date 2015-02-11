#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "ts.h"

#define MAXSTR 100

int m;
int n;
matrix_t *R;
matrix_t *C;

/* Terminate program and print error message from string 'info'  */
static void not_impl (char *info)
{
    fflush(stdout);
    fprintf(stderr, "Program terminated with error. %s\n", info);
    exit(1);
}

/* 
 * Reads a specified pattern 's' from stdin and terminates if the input
 * doesn't match the pattern.
 */
static void read_conf (const char *s)
{
    int i, len = strlen(s);
    char a[MAXSTR];
    for (i=0; i<len; i++)
    {
        if (scanf("%1s", &a[i])==EOF)
            /* make sure something is in  a  for error report       */
            a[i] = '\0';
        if (a[i] != s[i])
            /* the chars in  a  from stdin do not match those in s  */
        {
            fprintf(stderr, "\"%s\"  required from input, found \"%s\"\n",
                s, a);
            not_impl("");
        }
    }
}

/* Allocates memory for the bimatrix */
static void create_matrices()
{
    R = matrix_alloc(m, n);
    C = matrix_alloc(m, n);
}

/* Frees memory used to store the bimatrix */
static void free_matrices()
{
    matrix_free(R);
    matrix_free(C);
}

/* 
 * Reads the bimatrix from file using the gambit file format. 
 * Code originally from Codenotti et. al. 2008. 
 */
static void read_bimatrix_from_file(FILE *f, int* rdim1, int* rdim2) 
{
    char *buf = (char *) malloc(100 * sizeof(char));
    char c;
    int i, j, tmpn;
    size_t num_bytes = 100;
    int dim1, dim2; 
    double n1, n2;

    /* This checks if we are parsing the correct file type. */
    getline(&buf,&num_bytes,f);
    if( strncmp(buf,"NFG 1 D",(size_t)7) != 0 ) {
        fprintf(stderr,"NFG file corrupted, aborting\n");
        exit(1);
    }

    /* First we need to ignore all comments (and player names) from the .nfg file. */

    tmpn = 0;
    while(tmpn < 2) { /* This is to ignore comments */
        c = fgetc(f);
        if( c == '\"' )
            tmpn++;
    }
    tmpn = 0;
    while(tmpn < 2) { /* And this ignores player names */
        c = fgetc(f);
        if( c == '{' || c == '}')
            tmpn++;
    }

    fscanf(f,"%d %d",&dim1,&dim2);
    *rdim1 = dim1; *rdim2 = dim2;
    create_matrices();
    fgetc(f); fgetc(f);

    for(i = 0; i < dim2; i++) {
        for(j = 0; j < dim1; j++) {
            fscanf(f,"%lf %lf ",&n1,&n2);
            R->data[j][i] = n1; 
            C->data[j][i] =  n2;
        }
    }

    free(buf);
}

/* 
 * Read bimatrix game (R, C) of size i x j from stdin using the format
 * m= i
 * n= j
 * A= R
 * B= C
 */
static void read_bimatrix()
{
    read_conf("m=");
    scanf("%d", &m);
    read_conf("n=");
    scanf("%d", &n);

    create_matrices();

    read_conf("A=");
    matrix_read(R, stdin);
    read_conf("B=");
    matrix_read(C, stdin);
}

/* 
 * ./ts [i:d:ks:tw] 
 * Flags:
 * - '-i file' : Reads the bimatrix game from a file 'filename' using the 
 *   gambit file format otherwise reads in bimatrix from stdin.
 * - '-d delta': Sets the value of delta for TS to use. This should be 
 *   within the range (0, 1], otherwise default value of 0.1 is used.
 * - 'k' : UNUSED AT THE MOMENT. 
 * - 's seed' : Sets the random seed
 * - 't' : UNUSED AT THE MOMENT.
 */
int main(int argc, char **argv)
{
    int c;
    double delta = 0.1;

    char *nfg_file = NULL;

    while ((c = getopt(argc, argv, "i:wd:ks:t")) != -1)
        switch (c)
        {
            case 'i':
                nfg_file = optarg;
                break;
            case 'd':
                delta = atof(optarg);
                break;
            case 't':
                break;
            case 's':
                srand(atoi(optarg));
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

    if (delta <= 0 || delta > 1)
        delta = 0.1;

    if(nfg_file == NULL) {
        read_bimatrix();
    }
    else {
        FILE *f = fopen(nfg_file, "r");
        read_bimatrix_from_file(f, &m, &n);
        fclose(f);
    }

    run_ts(R, C, delta);

    free_matrices();
    return 0;
}
