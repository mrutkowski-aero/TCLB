#ifndef BXSPLINE_H
#define BXSPLINE_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
  extern "C" {
#endif

#define BXS_PERIODIC 0x01
#define BXS_NOEXTEND 0x02
    


double bxspline(double x, double s[], int order, double tab[], int Pars, bool per, int d);
void  vbxspline(double x, double s[], int order, int Pars, bool per, int d, double bxsrow[]);
double bxs_base(double s[], int order, double x, int i);
inline int mod(int i, int order, bool per);

#ifdef __cplusplus
  } // extern "C"
#endif

#endif // BOXSPLINE_H
