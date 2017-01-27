#include <stdio.h>
#include <string.h>
#include "implementation.h"

#include "exercise9.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

#define MAX_IT 30000
#define NEUMANN

void import_mesh_adcirc(mesh * m, char const * filename)
{
  importMeshAdcirc(m,filename);
}

void init_matrix(crs_matrix * mat, mesh const * m)
{
  initMatrix(mat,m);
}

void get_local_stiffness(double local_stiffness[3][3], mesh const * m,
                         size_t element_id)
{
  getLocalStiffness(local_stiffness,m, element_id);
}

void get_local_load(double local_load[3], mesh const * m, size_t element_id,
                    double (*fn_f)(double, double))
{
  getLocalLoad(local_load, m, element_id, fn_f);
}

void assemble_local2global_stiffness(double local_stiffness[3][3],
                                     crs_matrix * mat, mesh const * m,
                                     size_t element_id)
{
  assembleLocal2globalStiffness(local_stiffness, mat, m, element_id);
}

void assemble_local2global_load(double local_load[3], double * rhs,
                                mesh const * m, size_t element_id)
{
  assembleLocal2globalLoad(local_load, rhs,m,element_id);
}

void apply_dbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_g)(unsigned char, double, double))
{
  applyDbc(mat,rhs,m,fn_g);
}

void apply_nbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_v)(unsigned char, double, double, double const *))
{
  applyNbc(mat,rhs,m,fn_v);
}

void solve(crs_matrix const * mat, double * u, double const * rhs)
{
  Solve(mat,u,rhs);
}
