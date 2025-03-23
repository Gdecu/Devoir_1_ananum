#include "devoir_1.h"
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

    //k++;
    int count = 0;
    printf("k = %d\n", k);
    for (int i = 0; i < n * (k+1); i++) {
            /*if ((i) % (k+1) == 0) {
                A[i] = 0.0;
                count++;
            } else {
                A[i] = i - count;
            }*/
           A[i] = i;
    }
    
    //printf("Matrice bande initiale:\n");
    //print_band_matrix(A, n, k);

    //printf("Stockage de A :\n");
    //print_band_matrix_stockage(A, n, k);

    double *A_sym = (double *)malloc((n * (n + 1) / 2) * sizeof(double));
    symBand_to_sym(A, A_sym, n, k);

    printf("Matrice bande notation symmetrique:\n");
    print_sym_matrix(A_sym, n, k);

    tridiagonalize(A, n, k, d, e);
    
    //printf("Matrice après tridiagonalisation:\n");
    //print_band_matrix(A, n, k);
    
    printf("Diagonale principale:\n");
    for (int i = 0; i < n; i++) printf("%8.4f ", d[i]);
    printf("\n\nSous-diagonale:\n");
    for (int i = 0; i < n - 1; i++) printf("%8.4f ", e[i]);
    printf("\n");

    free(A);
    free(d);
    free(e);
}

int main( int argc, char *argv[] ) {

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);

    test_tridiagonalize(n, k);

    /* test pour qr eigs
    
    
    double *a = (double*)malloc(9 * sizeof(double));

    for (int i = 0; i < 9; i++) {
        a[i] = i;
    }

    printf("a_loc = %p, a[0] = %f, a[1] = %f, a[8] = %f\n", a, a[0], a[1], a[8]);

    a[0] = 10;
    a++;

    printf("a_loc = %p, a[0] = %f, a[1] = %f, a[8] = %f\n", a, a[0], a[1], a[8]);

    a[0] = 20;
    a--;
    a[0] = 30;

    printf("a_loc = %p, a[0] = %f, a[1] = %f, a[8] = %f\n", a, a[0], a[1], a[8]);

    */

    return 0;
}