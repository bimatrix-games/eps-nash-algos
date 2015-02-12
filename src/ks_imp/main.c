#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "ks.h"
#include "../utils/io.h"

/* 
 * ./ks+ [i:s:] 
 * Flags:
 * - '-i file' : Reads the bimatrix game from a file 'filename' using the 
 *   gambit file format otherwise reads in bimatrix from stdin.
 * - 's seed' : Sets the random seed
 */
int main(int argc, char **argv)
{
    int c, m, n;

    char *nfg_file = NULL;
    matrix_t **bimatrix;

    while ((c = getopt(argc, argv, "i:ws:")) != -1)
        switch (c)
        {
            case 'i':
                nfg_file = optarg;
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

    if(nfg_file == NULL) {
        bimatrix = read_bimatrix();
    }
    else {
        FILE *f = fopen(nfg_file, "r");
        bimatrix = read_bimatrix_from_file(f, &m, &n);
        fclose(f);
    }

    run_ks(bimatrix[0], bimatrix[1]);
    free_matrices(bimatrix);

    return 0;
}
