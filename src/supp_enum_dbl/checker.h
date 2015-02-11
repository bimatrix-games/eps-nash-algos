#include <stdio.h>
#include <stdlib.h>
#include "../utils/matrix.h"

#ifndef CHECKER_H
#define CHECKER_H

extern matrix_t *R;
extern matrix_t *C;
extern matrix_t *x;
extern matrix_t *y;
extern matrix_t *supp_x;
extern matrix_t *supp_y;

double best_response(matrix_t *mat, matrix_t *s);
double get_payoff(matrix_t *mat, matrix_t *x, matrix_t *y);
double compute_epsilon();
double compute_epsilon_supp();
int is_nash_support();

#endif
