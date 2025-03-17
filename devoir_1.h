#ifndef CODE2_H
#define CODE2_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


void tridiagonalize(double *A, int n, int k, double *d, double *e);
int step_qr_tridiag(double *d, double *e, double m, double eps);
int qr_eigs_(double *A, int n, int k, double eps, int max_iter, double *d);

#endif
