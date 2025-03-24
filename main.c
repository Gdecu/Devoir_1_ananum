#include "devoir_1.h"
#include "analyse.h"
 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <string.h> 
 
 #define idxSym(i, j) ((j) * ((j) + 1) / 2 + (i)) // symmetric
 #define idxSymBand(i, j, k) ((k) * (i+1) + j ) // band
 #define min_int(a, b) ((a) < (b) ? (a) : (b))
 
 void print_band_matrix(double *A, int n, int k) {
     for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
             if (fabs(i - j) <= k) {
                 if (i >= j) {
                     printf(" %8.4f ", A[idxSymBand(i, j, k)]);
                 } else {
                     printf(" %8.4f ", A[idxSymBand(j, i, k)]);
                 }
             } else {
                 printf(" %8.4f ", 0.0);
             }
         }
         printf("\n");
     }
     printf("\n");
 }
 
 void print_sym_matrix(double *A, int n, int k) {
     for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
             if (i >= j) {
                 printf(" %8.4f ", A[idxSym(j, i)]); // Accès à la partie triangulaire inférieure
             } else {
                 printf(" %8.4f ", A[idxSym(i, j)]); // Symétrie : A[i, j] = A[j, i]
             }
         }
         printf("\n");
     }
     printf("\n");
 }
 
 void print_band_matrix_stockage(double *A, int n, int k) {
     for (int i = 0; i < n * (k+1); i++) {
         printf("%8.4f ", A[i]);
         if ((i + 1) % (k + 1) == 0) {
             printf("\n");
         }
     }
     printf("\n");
 }
 
 void print_matrix_tridiagonal(double *d, double *e, int n) {
     for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
             if (i == j) {
                 printf(" %8.4f ", d[i]); 
             } else if (i == j - 1 || i == j + 1) {
                 printf(" %8.4f ", e[i < j ? i + 1 : j + 1]);  
             } else {
                 printf(" %8.4f ", 0.0); 
             }
         }
         printf("\n");
     }
     printf("\n");
 }
 
 void print_vector(double *v, int size) {
     for (int i = 0; i < size; i++) {
         printf("%lf ", v[i]);
     }
     printf("\n");
 }
 
 void symBand_to_sym( double *A_band, double *A_sym, int n, int k) {
     int index = 0;
     for (int i = 0; i < n; i++) {
         for (int j = 0; j <= i; j++) {  // On ne parcourt que la partie triangulaire inférieure
             if (j >= i - k) {  // Si l'élément appartient à la bande
                 A_sym[index++] = A_band[idxSymBand(i, j, k)];
             } else {
                 A_sym[index++] = 0.0;  // Optionnel : remplir avec 0 si hors de la bande
             }
         }
     }
 }
 
 void test_tridiagonalize(int n, int k) {
     double *A = (double *)malloc(n * (k + 1) * sizeof(double));
     double *d = (double *)malloc(n * sizeof(double));
     double *e = (double *)malloc((n - 1) * sizeof(double));
 
     int count = 0;
     printf("k = %d\n", k);
     for (int i = 0; i < n * (k+1); i++) {
            A[i] = i;
     }
     double *A_sym = (double *)malloc((n * (n + 1) / 2) * sizeof(double));
     symBand_to_sym(A, A_sym, n, k);
 
     printf("Matrice bande notation symmetrique:\n");
     print_sym_matrix(A_sym, n, k);
 
     tridiagonalize_band(A, n, k, d, e);
 
     printf("Matrice après tridiagonalisation:\n");
     print_matrix_tridiagonal(d, e, n);
 
     printf("Diagonale principale:\n");
     for (int i = 0; i < n; i++) printf("%8.4f ", d[i]);
     printf("\n\nSous-diagonale:\n");
     for (int i = 0; i < n - 1; i++) printf("%8.4f ", e[i]);
     printf("\n");
 
     free(A_sym);
     free(A);
     free(d);
     free(e);
 }
 
 void test_step_qr(int m, double eps) {
     double *d = (double *)malloc(m * sizeof(double));
     double *e = (double *)malloc(m * sizeof(double));
 
     if (d == NULL || e == NULL) {
         fprintf(stderr, "Memory allocation failed\n");
         exit(EXIT_FAILURE);
     }
 
     d[0] = 67.72;
     d[1] = 17.3154;
     d[2] = -70.5561;
     d[3] = 0.7163;
     e[0] = 0.7163;
     e[1] = -5.3277;
     e[2] = 289.1896;
     e[3] = 1.0;
     printf("Matrix before QR:\n");
     print_matrix_tridiagonal(d, e, m);
 
     step_qr_tridiag(d, e, m, eps);
 
     printf("Matrix after QR:\n");
     print_matrix_tridiagonal(d, e, m);
 
     free(d);
     free(e);
 }
 
 void test_qr_eigs_band(int n, int k, double eps, int max_iter){
     double *A = (double *)malloc(n * (k + 1) * sizeof(double));
     double *d = (double *)malloc(n * sizeof(double));
     if (A == NULL || d == NULL) {
         fprintf(stderr, "Memory allocation failed\n");
         exit(EXIT_FAILURE);
     }
     for (int i = 0; i < n; i++) {
         for (int j = 0; j <= k; j++) {
             A[i * (k + 1) + j] = ((double)rand() / RAND_MAX) * 10.0; 
         }
     }
     //printf("Matrice bande initiale:\n");
     //print_band_matrix(A, n, k);
     int a  = qr_eigs_band(A, n, k, eps, max_iter, d);
     //printf("Nombre d'itérations: %d\n", a);
     //printf("Valeur propre:\n");
     //print_vector(d, n);
     save_vector("eigenvalues.txt", d, n);
     free(A);
     free(d);
 }
 
 int main(int argc, char *argv[]) {
     if (argc != 3) {
         fprintf(stderr, "Usage: %s <m> <k>\n", argv[0]);
         return EXIT_FAILURE;
     }
 
     int n = atoi(argv[1]);
     int k = atoi(argv[2]);
 
     if (n < 2) {
         fprintf(stderr, "m must be at least 2\n");
         return EXIT_FAILURE;
     }
 
     //test_step_qr(n, 1e-6);
     //test_tridiagonalize(n, k);
     test_qr_eigs_band(n, k, 1e-6, 1000);


     double lx = 4.0;
     double ly = 2.0;
     int nx = 2;
     int ny = 2;
     n = nx * ny;
 
     double *A = create_matrix(nx, ny, lx, ly);
     double *d = (double *)malloc(n * sizeof(double));

    int a  = qr_eigs_band(A, n, nx, 1e-6, 1000, d);
    printf("Nombre d'itérations: %d\n", a);
    save_vector("eigenvalues_Laplacian.txt", d, n);


    free(A);
    free(d);

 
     return EXIT_SUCCESS;
 }