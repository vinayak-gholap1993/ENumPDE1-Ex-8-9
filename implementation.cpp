#include <cstdio>
#include <cstring>
#include "implementation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "exercise9.h"
#include "ordered_list.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

#define MAX_IT 30000
#define NEUMANN



extern "C" void importMeshAdcirc(mesh * m, char const * filename)
{
    std::string line;
    std::stringstream streamline;
    std::ifstream inputfile(filename);

    size_t NE = 0, NP =0;

    std::getline(inputfile, line); // read name of the file
    streamline << line;
    streamline >> m->name;
    streamline.clear();
    streamline.str(" ");

    std::getline(inputfile,line);// read the data NE , NP
    streamline << line;
    streamline >> NE >> NP;
    streamline.clear();
    streamline.str(" ");

    m->coords = new double [2 * NP];
    m->t2v = new size_t [3 * NE];
    m->id_v = new unsigned char [NP](); // default initialisation
    m->n_vertices = NP;
    m->n_triangles = NE;

    size_t n = 0;
    double x = 0.0, y =0.0;

    //Setting the coord data

    for( size_t i = 0; i < NP; ++i)
    {
        streamline.str(" ");
        streamline.clear();
        std::getline(inputfile, line);
        streamline << line;
        streamline >> n >> x >> y;
        m->coords[2 * i] = x;
        m->coords[2 * i +1] = y;
    }

    size_t node1, node2, node3, elementType, currentEleNum;

    // Setting the element connectivity vertices

    for(size_t k =0; k < NE; ++k)
    {
        streamline.str(" ");
        streamline.clear();
        std::getline(inputfile, line);
        streamline << line;
        streamline >> currentEleNum >> elementType >> node1 >> node2 >> node3;

        m->t2v[3 * k] = node1 -1; // As stated in pdf one eleNum to be subtracted to fit C style numbering
        m->t2v[3 * k + 1] = node2 -1;
        m->t2v[3 * k + 2] = node3 -1;
    }

    // Setting NOPE

    size_t NOPE = 0; // number of elevation specified boundary forcing segments
    streamline.str(" ");
    streamline.clear();
    std::getline(inputfile, line);
    streamline << line;
    streamline >> NOPE;

    //Setting NETA

    size_t NETA = 0; // total number of elevation specified boundary nodes
    streamline.str(" ");
    streamline.clear();
    std::getline(inputfile, line);
    streamline << line;
    streamline >> NETA;

    // Neumann BC
    // Setting num of nodes on open boundary segment

    for(size_t i = 0; i < NOPE; ++i )
    {
        size_t num_open_boundary_seg_nodes = 0, vertex = 0;
        streamline.str(" ");
        streamline.clear();
        std::getline(inputfile, line);
        streamline << line;
        streamline >> num_open_boundary_seg_nodes;

        //Setting vertex of Neumann BC

        for(size_t j = 0; j < num_open_boundary_seg_nodes; ++j)
        {
            streamline.str(" ");
            streamline.clear();
            std::getline(inputfile, line);
            streamline << line;
            streamline >> vertex; // Vertex connectivity of NBC

            m->id_v[vertex - 1] = i +1; // Again here Vertex-1 to be with C style numbering
        }
    }

    //Setting number of DBC segments

    size_t NBOU = 0;
    streamline.str(" ");
    streamline.clear();
    std::getline(inputfile, line);
    streamline << line;
    streamline >> NBOU;

    size_t NVEL = 0; // Total num of land boundary nodes
    streamline.str(" ");
    streamline.clear();
    std::getline(inputfile, line);
    streamline << line;
    streamline >> NVEL;

    size_t seg_index = 0, vertex =0, seg_nodes = 0;

    for(size_t i =0; i < NBOU; ++i)
    {
        streamline.str(" ");
        streamline.clear();
        std::getline(inputfile, line);
        streamline << line;
        streamline >> seg_nodes >> seg_index; // The num of Land or Island Segment nodes with the index for the segment

        for(size_t j =0; j < seg_nodes; ++j)
        {
            streamline.str(" ");
            streamline.clear();
            std::getline(inputfile, line);
            streamline << line;
            streamline >> vertex; // Vertex connectivity of DBC for land or island segment

            m->id_v[vertex - 1 ] = 100 + seg_index; // As given in PDF
        }
    }

    streamline.clear();
    inputfile.close();

}

extern "C" void getMesh(mesh * m, size_t n)
{
  /* Call to C++-routine */
  /* getMesh(m, n); */

   size_t coord_index, triangle_index;
   double x,y;
   size_t i,j;

   double h = 1.0/ n;
   double hsquare = h * h;
   double h_i[2];
   h_i[0] = h + hsquare;
   h_i[1] = h - hsquare;

   m->n_vertices = (n+1) * (n+1);
   m->n_triangles = 2 * n * n;

   //memory allocation of struct members

   m->coords = new double [2 * m->n_vertices];
   m->t2v = new size_t [3 * m->n_triangles];
   m->id_v = new unsigned char [m->n_vertices];

    for(i =0, y = 0.0; i <= n; ++i)
    {
        for(j=0, x = 0.0; j<= n ; ++j)
        {
            coord_index = i * (n+1) + j;        //lexicographic ordering of nodal points

            m->coords[2 * coord_index] = x;
            m->coords[2 * coord_index +1] = y;

            m->id_v[coord_index] = (j == 0) ? 4 : (j == n) ? 2 : (i == 0) ? 1 : (i == n) ? 3 : 0;

            x += h_i[j % 2];
        }
        y += h_i[i % 2];
    }

    for(i=0; i < n; ++i)
    {
        for(j=0; j < n; ++j )
        {
            triangle_index = 2 * (i * n+ j);

            // each quadrilateral will have 2 triangles
            // lower right triangle and upper left triangle

            //upper left triangle filling for t2v with coords (i,j), (i+1,j+1), (i+1,j) with counter clockwise ordering

            m->t2v[3 * triangle_index] = i * (n+1) + j;
            m->t2v[3 * triangle_index + 1] = (i+1) * (n+1) + (j+1);
            m->t2v[3 * triangle_index + 2] = (i+1) * (n+1) + j;

            //lower right triangle filling for t2v with coords (i,j), (i,j+1), (i+1,j+1) with counter clockwise ordering

            ++triangle_index;

            m->t2v[3 * triangle_index] = i * (n+1) + j;
            m->t2v[3 * triangle_index + 1] = i * (n+1) + (j+1);
            m->t2v[3 * triangle_index + 2] = (i+1) * (n+1) + (j+1);


        }
    }

}

extern "C" void initMatrix(crs_matrix * mat, mesh const * m)
{

    size_t vertex, triangle;
    size_t nnz = 0, ind;
    ordered_list * vertex2vertex_lists;
    ordered_list_element * iterator;

    vertex2vertex_lists = new ordered_list [m->n_vertices];
    if (vertex2vertex_lists == NULL)
        err_exit("Lists allocation failed");

    for (vertex = 0; vertex < m->n_vertices; ++vertex)
    {
      init_list(&vertex2vertex_lists[vertex]);
    }

    for (triangle = 0; triangle < m->n_triangles; ++triangle)
    {

      for (size_t i = 0; i < 3; ++i)
      {
        vertex = m->t2v[3 * triangle + i];

        for (size_t j = 0; j < 3; ++j)
        {
          insert_list(&vertex2vertex_lists[vertex], m->t2v[3 * triangle + j]);
        }
      }
    }

    for (vertex = 0; vertex < m->n_vertices; ++vertex)
    {
      nnz += vertex2vertex_lists[vertex].length;
    }

    mat->n_rows = mat->n_cols = m->n_vertices;
    mat->val = new double [nnz];
    mat->rowPtr = new size_t  [mat->n_rows + 1];
    mat->colInd = new size_t  [nnz];
    if (mat->val == NULL || mat->rowPtr == NULL || mat->colInd == NULL)
      err_exit("CRS-matrix allocation failed");

    //memset(mat->val, 0, sizeof(double) * nnz);
    mat->val = new double [nnz]();

    mat->rowPtr[0] = 0;

      for (vertex = 0; vertex < m->n_vertices; ++vertex)
      {
        ind = mat->rowPtr[vertex];
        iterator = vertex2vertex_lists[vertex].begin;
        while (iterator != NULL)
        {
          mat->colInd[ind++] = iterator->val;
          iterator = iterator->next;
        }
        mat->rowPtr[vertex+1] = ind;
      }
    //delete vertex2vertex_lists;
}

extern "C" void getLocalStiffness(double local_stiffness[3][3], mesh const * m,
                         size_t element_id)
{

    static const double S_1[3][3] = { {  0.5 , -0.5 , 0.0 },{ -0.5 ,  0.5 , 0.0 },{  0.0 ,  0.0 , 0.0 } };
    static const double S_2[3][3] = { {  1.0 , -0.5 , -0.5 }, { -0.5 ,  0.0 ,  0.5 }, { -0.5 ,  0.5 ,  0.0 } };
    static const double S_3[3][3] = { {  0.5 ,  0.0 , -0.5 }, {  0.0 ,  0.0 ,  0.0 }, { -0.5 ,  0.0 ,  0.5 } };

    double g1 = 0.0, g2 = 0.0, g3 = 0.0;
    double det_B_inv = 0.0;
    double * local_verts[3] ;
    size_t i, j;

    for(i = 0; i < 3; ++i)
    {
      local_verts[i] = &m->coords[2 * m->t2v[3 * element_id + i]];

      //std::cout<<"ele id: "<<element_id<<std::endl;
      //std::cout<<"t2v: "<<m->t2v[3 * element_id + i]<<std::endl;
      //std::cout<<"m coord: "<<&m->coords[2 * m->t2v[3 * element_id + i]]<<std::endl;


      //std::cout<<"Local verts at i: "<<local_verts[i][0]<<std::endl;
      //std::cout<<"Local verts at i: "<<local_verts[i][1]<<std::endl;

    }

    det_B_inv = 1.0 / fabs( (local_verts[1][0] - local_verts[0][0]) *(local_verts[2][1] - local_verts[0][1]) -
                           (local_verts[1][1] - local_verts[0][1]) *(local_verts[2][0] - local_verts[0][0]) );


    //std::cout<<"Det B inverse: "<<det_B_inv<<std::endl;

    g1 =  det_B_inv * ( (local_verts[2][0] - local_verts[0][0]) *(local_verts[2][0] - local_verts[0][0]) +
                             (local_verts[2][1] - local_verts[0][1]) *(local_verts[2][1] - local_verts[0][1]) );
    g2 = -det_B_inv * ( (local_verts[1][0] - local_verts[0][0]) * (local_verts[2][0] - local_verts[0][0]) +
                             (local_verts[1][1] - local_verts[0][1]) *(local_verts[2][1] - local_verts[0][1]) );
    g3 =  det_B_inv * ( (local_verts[1][0] - local_verts[0][0]) * (local_verts[1][0] - local_verts[0][0]) +
                             (local_verts[1][1] - local_verts[0][1]) *(local_verts[1][1] - local_verts[0][1]) );

    //std::cout<<"g1: "<<g1<<" "<<"g2: "<<g2<<" "<<"g3: "<<g3<<std::endl;

    for (i = 0; i < 3; ++i)
    {
      for (j = 0; j < 3; ++j)
      {
        local_stiffness[i][j] = g1 * S_1[i][j] + g2 * S_2[i][j] + g3 * S_3[i][j];
      }
    }
}

extern "C" void getLocalLoad(double local_load[3], mesh const * m, size_t element_id,
                    double (*fn_f)(double, double))
{

  double det_B;
  double * local_verts[3];
  size_t i;

  for(i = 0; i < 3; ++i)
  {
    local_verts[i] = &m->coords[2 * m->t2v[3 * element_id + i]];
  }

  det_B = fabs( (local_verts[1][0] - local_verts[0][0]) *(local_verts[2][1] - local_verts[0][1]) -
                (local_verts[1][1] - local_verts[0][1]) *(local_verts[2][0] - local_verts[0][0]) );

  // Applying Trapezoidal rule
  for (i = 0; i < 3; ++i)
  {
    local_load[i] = det_B / 6.0 * fn_f(local_verts[i][0], local_verts[i][1]);
  }
}

extern "C" void assembleLocal2globalStiffness(double local_stiffness[3][3],
                                     crs_matrix * mat, mesh const * m,
                                     size_t element_id)
{

  size_t i, j, ind;
  size_t global_i, global_j;

  for (i = 0; i < 3; ++i)
  {
    global_i = m->t2v[3 * element_id + i];

    for (j = 0; j < 3; ++j)
    {
      global_j = m->t2v[3 * element_id + j];

      for (ind = mat->rowPtr[global_i]; ind < mat->rowPtr[global_i + 1]; ++ind)
      {
        if (mat->colInd[ind] == global_j)
          break;
      }

      mat->val[ind] += local_stiffness[i][j];
    }
  }
}

extern "C" void assembleLocal2globalLoad(double local_load[3], double * rhs,
                                mesh const * m, size_t element_id)
{

  size_t i, global_i;

  if (rhs == NULL)
      err_exit("Rhs not given");

  for (i = 0; i < 3; ++i)
  {
    global_i = m->t2v[3 * element_id + i];

    rhs[global_i] += local_load[i];
  }
}

extern "C" void applyDbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_g)(unsigned char, double, double))
{
    //loop over vertices
    size_t j=0,col=0;
    double x=0,y=0;
    //std::cout<<m->n_vertices<<std::endl;

    for (size_t row=0; row<m->n_vertices; row++)
    {
            //find x,y coord of the vertex
            x=m->coords[2*row]; y=m->coords[2*row+1];

            //identify boundary nodes and check if it is not inside the domain

            //Should change the 0 to >= 100 in case of NEUMANN to avoid DBC

            if (m->id_v[row] > 0 )
            {
                    //edit matrix with 1.0 only on diagonal for DBC

                    for (j=mat->rowPtr[row]; j<mat->rowPtr[row+1]; j++)
                    {
                            col = mat->colInd[j];
                            //if colIndex=row number, data value = 1.0 otherwise 0.0
                            if (col == row)
                            {
                                    mat->val[j]=1.0;
                            }
                            else
                            {
                                    mat->val[j]=0.0;
                            }

                    }

                    //edit force vector with dirichlet value fn_g(row,x,y)
                    rhs[row]=fn_g(m->id_v[row],x,y);

            }
    }
}


extern "C" void applyNbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_v)(unsigned char, double, double, double const *))
{

    size_t v0 = 0, v1 = 0, v2 =0;


    double distance = 0.0, dx = 0.0, dy = 0.0;
    double normal[2];

    // Going through all elements
    for(size_t i =0; i < m->n_triangles; ++i)
    {
        //Vertices of the triangle

        v0 = m->t2v[3 * i];
        v1 = m->t2v[3 * i + 1];
        v2 = m->t2v[3 * i + 2];

        mat->n_cols = mat->n_cols; // To avoid unused variable error

        // Combinations of edges of triangle need to be checked for NBC
        // edges with 0, 1; 1, 2; 2,0 are to be checked


        // Check if one vertex is NB and other is not interior point but any boundary point

        // is v0 NB and v1 any boundary ?
        if(m->id_v[v0] == 1 && m->id_v[v1] != 0)
        {
            dx = (m->coords[2 * v0] - m->coords[2 * v1]);
            dy = (m->coords[2 * v0 +1 ] - m->coords[2 * v1 + 1]);
            distance = sqrt(dx * dx + dy * dy);

            normal[0] = -dy / distance;
            normal[1] = dx / distance;

            // Trapezoidal rule to integrate the rhs function
            rhs[v0] += 0.5 * distance * fn_v(m->id_v[v0], m->coords[2* v0], m->coords[2* v0 + 1], normal);

        }
        // is v1 NB and v0 any boundary?
        if(m->id_v[v1] == 1 && m->id_v[v0] != 0)
        {
            dx = (m->coords[2 * v1] - m->coords[2 * v0]);
            dy = (m->coords[2 * v1 +1 ] - m->coords[2 * v0 + 1]);
            distance = sqrt(dx * dx + dy * dy);

            normal[0] = -dy / distance;
            normal[1] = dx / distance;

            // Trapezoidal rule to integrate the rhs function
            rhs[v1] += 0.5 * distance * fn_v(m->id_v[v1], m->coords[2* v1], m->coords[2* v1 + 1], normal);

        }

        // is v1 NB and v2 any boundary?
        if(m->id_v[v1] == 1 && m->id_v[v2] != 0)
        {
            dx = (m->coords[2 * v1] - m->coords[2 * v2]);
            dy = (m->coords[2 * v1 +1 ] - m->coords[2 * v2 + 1]);
            distance = sqrt(dx * dx + dy * dy);

            normal[0] = -dy / distance;
            normal[1] = dx / distance;

            // Trapezoidal rule to integrate the rhs function
            rhs[v1] += 0.5 * distance * fn_v(m->id_v[v1], m->coords[2* v1], m->coords[2* v1 + 1], normal);

        }

        // is v2 NB and v1 any boundary?
        if(m->id_v[v2] == 1 && m->id_v[v1] != 0)
        {
            dx = (m->coords[2 * v2] - m->coords[2 * v1]);
            dy = (m->coords[2 * v2 +1 ] - m->coords[2 * v1 + 1]);
            distance = sqrt(dx * dx + dy * dy);

            normal[0] = -dy / distance;
            normal[1] = dx / distance;

            // Trapezoidal rule to integrate the rhs function
            rhs[v2] += 0.5 * distance * fn_v(m->id_v[v2], m->coords[2* v2], m->coords[2* v2 + 1], normal);

        }

        // is v2 NB and v0 any boundary?
        if(m->id_v[v2] == 1 && m->id_v[v0] != 0)
        {
            dx = (m->coords[2 * v2] - m->coords[2 * v0]);
            dy = (m->coords[2 * v2 +1 ] - m->coords[2 * v0 + 1]);
            distance = sqrt(dx * dx + dy * dy);

            normal[0] = -dy / distance;
            normal[1] = dx / distance;

            // Trapezoidal rule to integrate the rhs function
            rhs[v2] += 0.5 * distance * fn_v(m->id_v[v2], m->coords[2* v2], m->coords[2* v2 + 1], normal);

        }

        // is v0 NB and v2 any boundary?
        if(m->id_v[v0] == 1 && m->id_v[v2] != 0)
        {
            dx = (m->coords[2 * v0] - m->coords[2 * v2]);
            dy = (m->coords[2 * v0 +1 ] - m->coords[2 * v2 + 1]);
            distance = sqrt(dx * dx + dy * dy);

            normal[0] = -dy / distance;
            normal[1] = dx / distance;

            // Trapezoidal rule to integrate the rhs function
            rhs[v0] += 0.5 * distance * fn_v(m->id_v[v0], m->coords[2* v0], m->coords[2* v0 + 1], normal);

        }
    }

}

extern "C" void Solve(crs_matrix const * mat, double * u, double const * rhs)
{

    size_t it, k, i, ind;
    double a_kk = 0.0;

    for (it = 0; it < MAX_IT; ++it)
    {
      for (k = 0; k < mat->n_rows; ++k)
      {
        u[k] = rhs[k];
        for (ind = mat->rowPtr[k]; ind < mat->rowPtr[k+1]; ++ind)
        {
          i = mat->colInd[ind];
          if (i == k)
          {
            a_kk = mat->val[ind];
          }
          else
          {
            u[k] -= mat->val[ind] * u[i];
          }
        }
        u[k] /= a_kk;
      }
    }
}
