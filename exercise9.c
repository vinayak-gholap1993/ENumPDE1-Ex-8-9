#include <stdio.h>
#include <math.h>

#include "exercise9.h"
#include "fem.h"

/* Choose the testcase:
 *  1: u = x + y (zero RHS, should be solved up to machine precision)
 *  2: u = sin(...)*sin(...) (non-zero DBC and RHS)
 */
#define TESTCASE 2

double u(double x, double y)
{
#if TESTCASE == 1
  return x + y;
#elif TESTCASE == 2
  return sin(M_PI * (x-3.67e3)/3.5e4) * sin(M_PI * (y-5.33e3)/5.e4);
#else
#error("Invalid value for TESTCASE!");
#endif
}

double f(double x, double y)
{
#if TESTCASE == 1
  return 0*(x + y);
#elif TESTCASE == 2
  return M_PI * M_PI * (1/1.225e9 + 1/2.5e9) * u(x,y);
#else
#error("Invalid value for TESTCASE!");
#endif
}

double g(unsigned char i, double x, double y)
{
  return u(x,y) + 0*i;
}

double v(unsigned char i, double x, double y, double const * normal)
{
#if TESTCASE == 1
  return normal[0] + normal[1] + 0*(x+y+i);
#elif TESTCASE == 2
  return 0 * i + M_PI * (
      cos(M_PI*(x-3.67e3)/3.5e4)*sin(M_PI*(y-5.33e3)/5.e4)/3.5e4 * normal[0] +
      sin(M_PI*(x-3.67e3)/3.5e4)*cos(M_PI*(y-5.33e3)/5.e4)/5.0e4 * normal[1] );
#else
#error("Invalid value for TESTCASE!");
#endif
}

int main()
{
  double errors[2];

  fem(errors, &f, &g, &v, &u);

  printf("+------------+------------+\n");
  printf("| l2_norm    | inf_norm   |\n");
  printf("+------------+------------+\n");
  printf("| %8.4E | %8.4E |\n", errors[0], errors[1]);
  printf("+------------+------------+\n");

  return 0;
}
