#include "devoir_1.h"
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 

void print_band_matrix(double *A, int n, int k) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i - k && j <= i) {
                printf("%8.4f ", A[idxBand(i, j, k)]);
            } else {
                printf("%8.4f ", 0.0);
            }
        }
        printf("\n");
    }
    printf("\n");
}

void test_tridiagonalize_band() {
    int n = 4;
    int k = 1;
    double A[4 * 2] = {
        4.0, 1.0,
        1.0, 3.0,
        2.0, 5.0,
        0.0, 2.0
    };
    double d[4], e[4] = {0};
    
    printf("Matrice bande initiale:\n");
    print_band_matrix(A, n, k);
    
    tridiagonalize_band(A, n, k, d, e);
    
    printf("Matrice après tridiagonalisation:\n");
    print_band_matrix(A, n, k);
    
    printf("Diagonale principale:\n");
    for (int i = 0; i < n; i++) printf("%8.4f ", d[i]);
    printf("\n\nSous-diagonale:\n");
    for (int i = 0; i < n - 1; i++) printf("%8.4f ", e[i]);
    printf("\n");
}

int main( int argc, char *argv[] ) {

    test_tridiagonalize_band();

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