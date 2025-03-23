
#include "devoir_1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <complex.h>


#define idxBand(i, j, k) ((k) * (i+1) + j ) // band
#define min_int(a, b) ((a) < (b) ? (a) : (b))
#define max_int(a, b) ((a) > (b) ? (a) : (b))
#define square(x) ((x) * (x))
#define max_double(a, b) ((a) > (b) ? (a) : (b))
#define nearest(a, b, c)  ((fabs((a) - (c)) < fabs((b) - (c))) ? (a) : (b))


/**
 * @brief tridiagonalise une matrice symétrique bande par transformations de similitudes
 * @param A entrées de la matrice dans un tableau Row-Major :
 *          - de taille n × n si [format] est full
 *          - de taille n × (k + 1) si [format] est band
 * @param n taille de la matrice
 * @param k nombre de sous-diagonales non-nulles
 * @param d tableau de taille n qui contient en sortie la diagonale principale de la matrice tridiagonale (Hessenberg symétrique)
 * @param e est un tableau de taille n qui contient en sortie la sous-diagonale de la matrice tridiagonale dans ses n − 1 premiers éléments
 */
void tridiagonalize(double *A, int n, int k, double *d, double *e) {
    double c, s, a, b, r;
    double temp;    
    // Rotation de Givens
    // A --> G A G* = H
    //printf("A --> A G\n");
    for (int j = 0; j < n-1; j++){
        for (int i = k + j - 1; i >= j + 1; i--){

            a = A[idxBand(i, j, k)];         // Élément diagonal
            b = A[idxBand(i+1, j, k)];       // Élément sous-diagonal
            //printf("a = %f, b = %f\n", a, b);
            //printf("i = %d, j = %d\n", i, j);

            if (b == 0) {continue;}             // Si l'élément sous-diagonal est nul, on passe à l'itération suivante

            r = sqrt(a * a + b * b);
            c = a / r;
            s = -b / r;
            //printf("c = %f, s = %f, r = %f\n", c, s, r);

            // Mise à jour des éléments de la matrice en appliquant la rotation

            double Aii = A[idxBand(i, i, k)];
            double Ai1i = A[idxBand(i+1, i, k)];
            double Ai1i1 = A[idxBand(i+1, i+1, k)];


            double Aij, Aij1; 
            int col = j + 1, row = i; int reverse = 0;
            // On mets à jour le reste
            for (int x = 0; x < k; x++){
                if (col == row) { // On ne touche pas aux éléments qui seront mis à jour plus bas
                    reverse = 1;
                    row += 2;
                    continue;
                } 
                if (row > n - 1) break;
                Aij = A[idxBand(row, col, k)]; //printf("Aij = %f\n", Aij);
                if (reverse == 0) {
                    Aij1 = A[idxBand(row + 1, col, k)];
                } else {
                    Aij1 = A[idxBand(row, col + 1, k)];
                } //printf("Aij1 = %f\n", Aij1);

                A[idxBand(row, col, k)] = c * Aij - s * Aij1; //printf("A[%d, %d], %f * %f - %f * %f = %f\n", row, col,c, Aij, s, Aij1, c * Aij - s * Aij1);
                A[idxBand(row, col + 1, k)] = s * Aij + c * Aij1; //printf("A[%d, %d], %f * %f + %f * %f = %f\n", row, col + 1, s, Aij, c, Aij1, s * Aij + c * Aij1);

                reverse == 0 ?  col++ : row++;
            }
            //temp = - s * A[idxBand(row, col, k)];  // élément non stockable dans la matrice symmetrique bande
            //A[idxBand(row, col, k)] = c * A[idxBand(row, col, k)];

            A[idxBand(i, j, k)] = r;                // A[i, j] = r
            A[idxBand(i+1, j, k)] = 0;              // A[i+1, j] = 0
            A[idxBand(i, i, k)] = c * c * Aii - 2 *c * s * Ai1i + s * s * Ai1i1;
            A[idxBand(i+1, i, k)] = c * s * Aii + (c*c - s*s) * Ai1i - c * s * Ai1i1;
            A[idxBand(i+1, i+1, k)] = s * s * Aii + 2 * c * s * Ai1i + c * c * Ai1i1;

        }
    }

    // On récupère les éléments diagonaux et sous-diagonaux de la matrice tridiagonale
    for (int i = 0; i < n; i++) {
        d[i] = A[idxBand(i, i, k)];
        if (i < n - 1) {
            e[i] = A[idxBand(i + 1, i, k)];
        }
    }
}

/**
 * @brief effectue une étape de l’algorithme QR avec un shift de Wilkinson µ sur une matrice tridiagonale symétrique
 * @param d tableau de taille n qui contient la diagonale principale de la matrice tridiagonale, mis à jour en sortie de telle sorte que H_{k+1} ← Q^{∗}_{k} H_{k} Q_{k}
 * @param e tableau de taille n qui contient la sous-diagonale de la matrice tridiagonale, mis à jour en sortie de telle sorte que H_{k+1} ← Q^{∗}_{k} H_{k} Q_{k}
 * @param m indique la taille de la matrice active (m ≤ n), c’est-à-dire sur laquelle on n’a pas encore isolé les valeurs propres
 * @param eps la tolérance pour déclarer qu’une valeur propre a été isolée : Golub et Van Loan proposent |hp+1,p| ≤ ϵ(|hpp| + |hp+1,p+1|).
 * @return nombre la nouvelle taille de la matrice active
 */
int step_qr_tridiag(double *d, double *e, int m, double eps) {
    if (m < 2) return 0;
    double mu = nearest((d[m-1] + d[m-2] - sqrt(square(d[m-1] + d[m-2]) + 4 * square(e[m-1])))/2, (d[m-1] + d[m-2] + sqrt(square(d[m-1] + d[m-2]) + 4 * square(e[m-1])))/2, d[m-1]);
    for (int i = 0; i < m - 1; i++) d[i] -= mu;
    double a,b,c,s,r,t_1,t_2;
    double cosinus[m-1];
    double sinus[m-1];
    t_2 = e[1];
    for (int i = 0; i <  m-1; i ++){
        a = d[i];
        b = t_2;
        r = hypot(a, b);
        c = a / r;
        s = -b / r;
        cosinus[i] = c;
        sinus[i] = s;
        if (i < m-2) t_2 = e[i + 2];
        t_1 = d[i + 1];
        d[i + 1] = s * e[i + 1] + c * d[i + 1];
        e[i + 1] = c * e[i + 1] - s * t_1;
        d[i] = r;
        if (i < m-2) e[i + 2] = c * e[i + 2];
    }
    for (int i = 0; i < m - 1; i++){
        c = cosinus[i];
        s = sinus[i];
        d[i] = c * d[i] - s * e[i+1];
        e[i+1] = - s * d[i+1];
        d[i+1] = c * d[i+1];
    }
    for (int i = 0; i < m - 1; i++) d[i] += mu;
    if (fabs(e[m-1]) <= eps * (fabs(d[m-2]) + fabs(d[m-1]))) {
        return m - 1;  
    }
    return m;
}

/** 
 * @brief calcule l’entièreté du spectre d’une matrice bande symétrique en faisant appel aux deux fonctions précedentes
 * @param A entrées de la matrice dans un tableau Row-Major :
 *          - de taille n × n si [format] est full
 *          - de taille n × (k + 1) si [format] est band
 * @param n taille de la matrice
 * @param k nombre de sous-diagonales non-nulles
 * @param eps la tolérance pour déclarer qu’une valeur propre a été isolée : Golub
 * @param max_iter nombre maximal d’itérations qr qu'on souhaite effectuer
 * @param d tableau de taille n qui contient en sortie les valeurs propres de la matrice A
 * @return retourne le nombre d’itérations nécessaires, ou bien -1 si l’algorithme n’a pas convergé
*/
int qr_eigs_(double *A, int n, int k, double eps, int max_iter, double *d) {
    double *e = (double *)malloc(n * sizeof(double));
    if (e == NULL) {
        free(e);
        return -1;
    }
    tridiagonalize(A, n, k, d, e);
    for (int i = 0; i < max_iter; i++) {
        int m = step_qr_tridiag(d, e, n, eps);
        if (m == 0) {
            free(e);
            return i;
        }
        n = m;
    }
    return -1;
}