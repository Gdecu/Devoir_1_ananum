#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "devoir_1.h"
#include <time.h>


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

    k = nx;
    alpha = 1. / dx2;
    beta = 1. / dy2;
    gamma = 2 * (alpha + beta);
    lda = k + 1;
    double *L = (double *)calloc(size * lda, sizeof(double));
    if (L == NULL) {
        printf("Error allocating memory for the matrix\n");
        exit(1);
    }
    for (int l = 0; l < size; l++) {
        L[idxSymBand(l,l,k)] = -beta; // (i,j)->(i-1,j)
        if (l % k != 0)
            L[idxSymBand(l,l-1,k)] = -alpha; // (i,j)->(i,j-1)
        L[idxSymBand(l,l-k,k)] = +gamma;     // (i,j)->(i  ,j)
    }    
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

void save_vector(char *name, double *d, int n) {
    FILE *f = fopen(name, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", name);
        exit(1);
    }
    fprintf(f, "# n = %d\n", n);
    for (int i = 0; i < n; i++) {
        fprintf(f, "%20.15le\n", d[i]);
    }
    fclose(f);
}

void time_complexity_qr_eig(double *L, double *d, double *times, int m, double lx, double ly) {
    int step = 1;
    int nb_repetitions = 10;

    int n, nx;
    double tmp;
    clock_t debut, fin;
    
    for (nx = 0; nx < m; nx += step) {
        n = nx * nx;
        L = create_matrix(nx, nx, lx, ly);

        tmp = 0.0;
        for (int i = 0; i < nb_repetitions; i++) {
            debut = clock();
            qr_eigs_band(L, n , nx, 1e-6, 500, d);
            fin = clock();
            tmp += ((double)(fin - debut));
        }
        double temps_moyen = tmp / (nb_repetitions * CLOCKS_PER_SEC);
        printf("n = %d, Temps moyen = %f secondes\n", n, temps_moyen);
        times[nx-1] = temps_moyen;        
    }
}