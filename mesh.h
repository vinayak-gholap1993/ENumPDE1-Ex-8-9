#ifndef MESH_H
#define MESH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

/*
 * A struct holding all relevant data structures to describe the triangulation.
 */
typedef struct {
  /* Interleaved node coordinates */
  double * coords;

  /* Vertex indices of each triangle in counter-clock-wise order */
  size_t * t2v;

  /* Vertex type id  (0: interior, other: boundaries) */
  unsigned char * id_v;

  /* Name of the mesh */
  char name[120];

  /* Number of vertices */
  size_t n_vertices;

  /* Number of triangles */
  size_t n_triangles;
} mesh;

/*
 * Generates a triangulation of the unit square for a given n.
 */
void get_mesh(mesh * m, size_t n);

/*
 * Imports a mesh in ADCIRC format from a text file.
 */
void import_mesh_adcirc(mesh * m, char const * filename);

/*
 * Checks whether the given mesh is allocated.
 */
unsigned char is_allocated_mesh(mesh const * m);

/*
 * Prints the content of the mesh data structures to stdout.
 */
void print_mesh(mesh const * m);

/*
 * Frees all allocated resources of a mesh (but not the mesh object itself!).
 */
void free_mesh(mesh * m);

#ifdef __cplusplus
}
#endif

#endif
