#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "bbm.h"
#include "../utils/io.h"

/* 
 * ./bbm [i:s:t] 
 * Flags:
 * - '-i file' : Reads the bimatrix game from a file 'filename' using the 
 *   gambit file format otherwise reads in bimatrix from stdin.
 *   within the range (0, 1], otherwise default value of 0.1 is used.
 * - 's seed' : Sets the random seed
 * - 't' : Use BBM2 instead of the default BBM1
 */

int main(int argc, char **argv)
{
    int c, m, n;
    char *nfg_file = NULL;
    int bbm = 0;
    matrix_t **bimatrix;

    while ((c = getopt(argc, argv, "i:wks:t")) != -1)
        switch (c)
        {
            case 'i':
                nfg_file = optarg;
                break;
            case 't':
                bbm = 1;
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

    run_bbm(bimatrix[0], bimatrix[1], bbm);
    free_matrices(bimatrix);

    return 0;
}
