#include "exercise9.h"

#ifdef __cplusplus
extern "C" void getMesh(mesh * m, size_t n);
#else
extern void getMesh(mesh * m, size_t n);
#endif

#ifdef __cplusplus
extern "C" void importMeshAdcirc(mesh * m, char const * filename);
#else
extern void importMeshAdcirc(mesh * m, char const * filename);
#endif

#ifdef __cplusplus
extern "C" void initMatrix(crs_matrix * mat, mesh const * m);
#else
extern void initMatrix(crs_matrix * mat, mesh const * m);
#endif

#ifdef __cplusplus
extern "C" void getLocalStiffness(double local_stiffness[3][3], mesh const * m,
                         size_t element_id);
#else
extern void getLocalStiffness(double local_stiffness[3][3], mesh const * m,
                         size_t element_id);
#endif

#ifdef __cplusplus
extern "C" void getLocalLoad(double local_load[3], mesh const * m, size_t element_id,
                    double (*fn_f)(double, double));
#else
extern void getLocalLoad(double local_load[3], mesh const * m, size_t element_id,
                    double (*fn_f)(double, double));
#endif

#ifdef __cplusplus
extern "C" void assembleLocal2globalStiffness(double local_stiffness[3][3],
                                     crs_matrix * mat, mesh const * m,
                                     size_t element_id);
#else
extern void assembleLocal2globalStiffness(double local_stiffness[3][3],
                                     crs_matrix * mat, mesh const * m,
                                     size_t element_id);
#endif

#ifdef __cplusplus
extern "C" void assembleLocal2globalLoad(double local_load[3], double * rhs,
                                mesh const * m, size_t element_id);
#else
extern void assembleLocal2globalLoad(double local_load[3], double * rhs,
                                mesh const * m, size_t element_id);
#endif

#ifdef __cplusplus
extern "C" void applyDbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_g)(unsigned char, double, double));
#else
extern void applyDbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_g)(unsigned char, double, double));
#endif

#ifdef __cplusplus
extern "C" void applyNbc(crs_matrix * mat, double * rhs, mesh const * m,
                          double (fn_v)(unsigned char, double, double, double const *));
#else
extern void applyNbc(crs_matrix * mat, double * rhs, mesh const * m,
                      double (fn_v)(unsigned char, double, double, double const *));
#endif

#ifdef __cplusplus
extern "C" void Solve(crs_matrix const * mat, double * u, double const * rhs);
#else
extern void Solve(crs_matrix const * mat, double * u, double const * rhs);
#endif
