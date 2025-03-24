#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


void tridiagonalize(double *A, int n, int k, double *d, double *e);
int step_qr_tridiag(double *d, double *e, int m, double eps);
int qr_eigs_(double *A, int n, int k, double eps, int max_iter, double *d);
void print_band_matrix(double *A, int n, int k);
void print_sym_matrix(double *A, int n, int k);
void symBand_to_sym( double *A_band, double *A_sym, int n, int k);
