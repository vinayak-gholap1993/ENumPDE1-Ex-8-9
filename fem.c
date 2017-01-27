#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fem.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

void fem(double errors[2], double (*fn_f)(double, double),
         double (*fn_g)(unsigned char, double, double),
         double (*fn_v)(unsigned char, double, double, double const *),
         double (*fn_u)(double, double))
{
  mesh m;
  crs_matrix mat;
  double * u, * rhs, * err, * u_exact;
  double local_stiffness[3][3];
  double local_load[3];
  size_t elem, i;

  /* 1. Import mesh */
  import_mesh_adcirc(&m, "fort.14");
  printf("Mesh: %s\nElements: %d\nVertices: %d\n",
      m.name, (int)m.n_triangles, (int)m.n_vertices);

#ifdef PRINT_DEBUG
  print_mesh(&m);
  write_vtk("mesh.vtk", &m, NULL);
#endif

  /* 2. Allocate the linear system */
  init_matrix(&mat, &m);

  rhs = (double *) malloc(sizeof(double) * m.n_vertices);
  if (rhs == NULL)
    err_exit("Allocation of right hand side failed!");
  memset(rhs, 0, sizeof(double) * m.n_vertices);

  /* 3. Assemble the matrix */
  for (elem = 0; elem <  m.n_triangles; ++elem)
  {
    /* Compute local stiffness and load */
    get_local_stiffness(local_stiffness, &m, elem);
    get_local_load(local_load, &m, elem, fn_f);

#ifdef PRINT_DEBUG
    print_local_stiffness(local_stiffness);
    print_local_load(local_load);
#endif

    /* insert into global matrix and rhs */
    assemble_local2global_stiffness(local_stiffness, &mat, &m, elem);
    assemble_local2global_load(local_load, rhs, &m, elem);
  }

#ifdef PRINT_DEBUG
  printf("Matrix after assembly:\n");
  print_matrix(&mat);
#endif

  /* 4. Apply boundary conditions */
  apply_nbc(&mat, rhs, &m, fn_v);
  apply_dbc(&mat, rhs, &m, fn_g);

#ifdef PRINT_DEBUG
  printf("Matrix after application of BCs:\n");
  print_matrix(&mat);
#endif

  /* 5. Solve the linear system */
  u = (double *) malloc(sizeof(double) * m.n_vertices);
  if (u == NULL) err_exit("Allocation of solution vector failed!");
  memset(u, 0, sizeof(double) * m.n_vertices);
  solve(&mat, u, rhs);

  /* 6. Evaluate analytical solution and pointwise error */
  u_exact = (double *) malloc(sizeof(double) * m.n_vertices);
  err = (double *) malloc(sizeof(double) * m.n_vertices);
  if (err == NULL || u_exact == NULL)
    err_exit("Allocation of error or exact solution vectors failed!");

  errors[0] = errors[1] = 0;
  for (i = 0; i < m.n_vertices; ++i)
  {
    u_exact[i] = fn_u(m.coords[2*i], m.coords[2*i+1]);
    err[i]     = u_exact[i] - u[i];

    errors[0] += err[i] * err[i];
    errors[1] = errors[1] > fabs(err[i]) ? errors[1] : fabs(err[i]);
  }

  /* 7. Visualize results */
  write_vtk("u.vtk", &m, u);
  write_vtk("u_exact.vtk", &m, u_exact);
  write_vtk("err.vtk", &m, err);

  /* free allocated resources */
  free(u);
  u = NULL;
  free(rhs);
  rhs = NULL;
  free(u_exact);
  u_exact = NULL;
  free(err);
  err = NULL;
  free_matrix(&mat);
  free_mesh(&m);
}

void write_vtk(char const * file_name, mesh const * m, double const * data)
{
  FILE * file;
  size_t i;

  if (!is_allocated_mesh(m)) err_exit("Mesh is not allocated!");

  /* Open file for writing */
  file = fopen(file_name, "w");
  if (file == NULL) err_exit("Failed to open file!");

  /* Write file header */
  fprintf(file, "# vtk DataFile Version 2.0\n");
  fprintf(file, "%s\n", m->name);
  fprintf(file, "ASCII\n");
  fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

  /* Write point coordinates */
  fprintf(file, "POINTS %d double\n", (int) m->n_vertices);
  for (i = 0; i < m->n_vertices; ++i)
  {
    fprintf(file, "%.3e %.3e %.3e\n",
        m->coords[2*i], m->coords[2*i + 1], 0.0);
  }

  /* Write cell connectivity */
  fprintf(file, "CELLS %d %d\n", (int) m->n_triangles, (int)(m->n_triangles*4));
  for (i = 0; i < m->n_triangles; ++i)
  {
    fprintf(file, "3 %d %d %d\n",
      (int) m->t2v[3*i], (int) m->t2v[3*i+1], (int) m->t2v[3*i+2]);
  }

  /* Write cell types */
  fprintf(file, "CELL_TYPES %d\n", (int) m->n_triangles);
  for (i = 0; i < m->n_triangles; ++i)
  {
    fprintf(file, "5\n");
  }

  /* Write point data */
  fprintf(file, "POINT_DATA %d\n", (int) m->n_vertices);

  /* data */
  if (data != NULL)
  {
    fprintf(file, "SCALARS data double 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (i = 0; i < m->n_vertices; ++i)
    {
      fprintf(file, "%.3e\n", data[i]);
    }
  }

  /* Vertex ids */
  fprintf(file, "SCALARS id_v int 1\n");
  fprintf(file, "LOOKUP_TABLE default\n");
  for (i = 0; i < m->n_vertices; ++i)
  {
    fprintf(file, "%d\n", (int) m->id_v[i]);
  }


  /* Close file */
  fclose(file);
}

void print_local_stiffness(double local_stiffness[3][3])
{
  size_t i, j;
  for (i = 0; i < 3; ++i)
  {
    for (j = 0; j < 3; ++j)
    {
      printf("%5.2f ", local_stiffness[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_local_load(double local_load[3])
{
  size_t i;
  for (i = 0; i < 3; ++i)
  {
    printf("%5.2f ", local_load[i]);
  }
  printf("\n\n");
}

