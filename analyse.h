#ifndef CODE2_H
#define CODE2_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double *create_matrix(int nx, int ny, double lx, double ly);
void print_sym_band(double *A, int n, int b, char *name);
void save_eigenvalues(char *name, double *d, int n);
#endif