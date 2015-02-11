#include <ilcplex/cplex.h>
#include "../utils/matrix.h"
#include "../utils/cplp.h"

typedef struct {
    matrix_t *x;
    matrix_t *y;
}sol_t;

void run_ks(matrix_t *R, matrix_t *C);
