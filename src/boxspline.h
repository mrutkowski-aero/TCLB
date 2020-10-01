#ifndef BOXSPLINE1_H
#define BOXSPLINE1_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
  extern "C" {
#endif

#define BS1_PERIODIC 0x01
#define BS1_NOEXTEND 0x02
    

double bs1_base(double s[], int k, double x, int i);
int mod(int i, int k);
double bs1_tab (double s[], int k, double val[], int n, double a, double b, double x, int d);
void bs1_design (double s[], int k, int n, double a, double b, double x, int d, double ret[]);

#ifdef __cplusplus
  } // extern "C"
#endif

#endif // BOXSPLINE1_H
