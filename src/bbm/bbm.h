#include <ilcplex/cplex.h>
#include "util.h"
#include "../utils/matrix.h"
#include "../utils/cplp.h"

typedef struct {
    matrix_t *x;
    matrix_t *y;
}bbm_sol_t;

void run_bbm(matrix_t *R, matrix_t *C, int bbm);
