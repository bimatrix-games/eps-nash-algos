#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "ts.h"
#include "../utils/io.h"

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
    int c, m, n;
    double delta = 0.1;

    char *nfg_file = NULL;
    matrix_t **bimatrix;

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
        bimatrix = read_bimatrix();
    }
    else {
        FILE *f = fopen(nfg_file, "r");
        bimatrix = read_bimatrix_from_file(f, &m, &n);
        fclose(f);
    }

    run_ts(bimatrix[0], bimatrix[1], delta);

    free_matrices(bimatrix);
    return 0;
}
