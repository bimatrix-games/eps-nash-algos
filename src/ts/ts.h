#include <ilcplex/cplex.h>
#include "../utils/util.h"
#include "../utils/matrix.h"
#include "../utils/cplp.h"

typedef struct {
    matrix_t *x;
    matrix_t *y;
}sol_t;

typedef struct {
    sol_t prob;
    cplp_sol_t *lp_sol;
}steep_sol_t;

typedef struct ts_sol_t{
    sol_t x_y;
    sol_t xp_yp;
    sol_t xh_yh;
    long iterations;
}ts_sol_t;

void run_ts(matrix_t *R, matrix_t *C, double delta);
