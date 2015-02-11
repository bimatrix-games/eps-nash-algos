#include <ilcplex/cplex.h>
#include "../utils/matrix.h"
#include "../utils/cplp.h"

typedef struct sol_t{
    matrix_t *x;
    matrix_t *y;
}sol_t;

typedef struct {
    matrix_t *R;
    matrix_t *C;
} ks_args;

void run_ks(matrix_t *R, matrix_t *C);
