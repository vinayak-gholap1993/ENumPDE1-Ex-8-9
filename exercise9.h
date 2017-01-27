#ifndef EXERCISE_FIVE_H
#define EXERCISE_FIVE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "crs_matrix.h"
#include "fem.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Evaluates the exact solution at point (x,y)
 */
double u(double x, double y);

/*
 * Evaluates the right hand side at point (x,y)
 */
double f(double x, double y);

/*
 * Returns Dirichlet boundary value g_i at point (x,y),
 */
double g(unsigned char i, double x, double y);

/*
 * Returns Neumann boundary value v_i at point (x,y)
 */
double v(unsigned char i, double x, double y, double const * normal);

#ifdef __cplusplus
}
#endif

#endif
