#include "io.h"
/* Terminate program and print error message from string 'info'  */
void not_impl (char *info)
{
    fflush(stdout);
    fprintf(stderr, "Program terminated with error. %s\n", info);
    exit(1);
}

/* 
 * Reads a specified pattern 's' from stdin and terminates if the input
 * doesn't match the pattern.
 */
void read_conf (const char *s)
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

/* 
 * Reads the bimatrix from file using the gambit file format. 
 * Code originally from Codenotti et. al. 2008. 
 */
matrix_t **read_bimatrix_from_file(FILE *f, int* rdim1, int* rdim2) 
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
    matrix_t **bimatrix = create_matrices();
    matrix_t *R = bimatrix[0];
    matrix_t *C = bimatrix[1];
    fgetc(f); fgetc(f);

    for(i = 0; i < dim2; i++) {
        for(j = 0; j < dim1; j++) {
            fscanf(f,"%lf %lf ",&n1,&n2);
            R->data[j][i] = n1; 
            C->data[j][i] =  n2;
        }
    }

    free(buf);
    return bimatrix;
}

/* 
 * Read bimatrix game (R, C) of size i x j from stdin using the format
 * m= i
 * n= j
 * A= R
 * B= C
 */
matrix_t **read_bimatrix()
{
    read_conf("m=");
    scanf("%d", &m);
    read_conf("n=");
    scanf("%d", &n);

    matrix_t **bimatrix = create_matrices();

    read_conf("A=");
    matrix_read(bimatrix[0], stdin);
    read_conf("B=");
    matrix_read(bimatrix[1], stdin);
    return bimatrix;
}

matrix_t **create_matrices(int m, int n)
{
    matrix_t **bimatrix = malloc(sizeof(matrix_t *) * 2);
    bimatrix[0] = matrix_alloc(m, n);
    bimatrix[1] = matrix_alloc(m, n);
    return bimatrix;
}

void free_matrices(matrix_t **bimatrix)
{
    matrix_free(bimatrix[0]);
    matrix_free(bimatrix[1]);
    free(bimatrix);
}
