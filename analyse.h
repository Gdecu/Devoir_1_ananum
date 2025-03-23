#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double *create_matrix( int nx, int ny, double lx, double ly);
void print_sym_band(double *A, int n, int b, char *name);
void save_vector(char *name, double *d, int n);
void time_complexity_qr_eig(double *L, double *d, double *times, int m, double lx, double ly);