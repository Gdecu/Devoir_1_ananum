#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "devoir_1.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SQUARE(a) ((a) * (a))
#define TOL_PRINT 0e-15
#define MAX_SIZE_PRINT 25
#define idxSymBand(i, j, k) ((k) * (i+1) + j ) // SYM band
#define idxSym(i, j) ((i) * ((i) + 1) / 2 + (j)) // symmetric#define sym_index(i, j) (j * (j + 1) / 2 + i)


/**
 * @brief Fills the matrix with the (opposite of the) Laplacian operator
 * @param nx Number of inner nodes in the x direction
 * @param ny Number of inner nodes in the y direction
 * @param lx Length of the domain in the x direction
 * @param ly Length of the domain in the y direction
 * @return the pointer to the discrete Laplacian matrix
 */
double *create_matrix(int nx, int ny, double lx, double ly) {

    int lda, k;
    int size = nx * ny; // Number of nodes/unknowns
    double dx2 = SQUARE(lx / (nx + 1));
    double dy2 = SQUARE(ly / (ny + 1));
    double alpha, beta, gamma;
    double *L;

    // Choice of node numbering, here nx=4, ny=5
    // j\i    0   1   2   3
    //     .  .   .   .   .  .
    // 0   .  0   1   2   3  .
    // 1   .  4   5   6   7  .
    // 2   .  8   9  10  11  .
    // 3   . 12  13  14  15  .
    // 4   . 16  17  18  19  .
    //     .  .   .   .   .  .

    k = nx;
    alpha = 1. / dx2;
    beta = 1. / dy2;
    gamma = 2 * (alpha + beta);
        lda = k + 1;
        L = (double *)calloc(size * lda, sizeof(double));
        for (int l = 0; l < size; l++) {
            L[l * lda + k - k] = -beta; // (i,j)->(i-1,j)
            if (l % k != 0)
                L[l * lda + k - 1] = -alpha; // (i,j)->(i,j-1)
            L[l * lda + k - 0] = +gamma;     // (i,j)->(i  ,j)
        }    
    return L;
}

void print_sym_band(double *A, int n, int b, char *name) {
    if (n > MAX_SIZE_PRINT)
        return;
    printf("\nSymmetric band matrix %s\n", name);
    int idx;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i - b; j++) {
            printf("%6s ", "");
        }
        for (int j = MAX(0, b - i); j < b + 1; j++) {
            idx = i * (b + 1) + j;
            // printf("%6d ", idx);
            if (fabs(A[idx]) < TOL_PRINT)
                printf("%6s ", "");
            else
                printf("%6.2lf ", A[idx]);
        }
        printf("\n");
    }
    printf("\n");
}

void save_eigenvalues(char *name, double *d, int n) {
    FILE *f = fopen(name, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", name);
        exit(1);
    }
    fprintf(f, "# Eigenvalues\n");
    fprintf(f, "# n = %d\n", n);
    for (int i = 0; i < n; i++) {
        fprintf(f, "%20.15le\n", d[i]);
    }
    fclose(f);
}


/**
 * @brief Computes eigenvalues and eigenvectors of a discrete Laplacian matrix
 * @param lx Length of the domain in the x direction
 * @param ly Length of the domain in the y direction
 * @param nx Number of inner nodes in the x direction
 * @param ny Number of inner nodes in the y direction
 */
