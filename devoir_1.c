#include "devoir_1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 

#define idxBand(i, j, k) ((i) * ((k) + 1) + ((j) - (i) + (k))) // band
#define min_int(a, b) ((a) < (b) ? (a) : (b))
#define max_int(a, b) ((a) > (b) ? (a) : (b))
#define square(x) ((x) * (x))
#define max_double(a, b) ((a) > (b) ? (a) : (b))

/*

 double square( double x){
    return x * x;
}
double max_double(double a, double b) {
    return (a > b) ? a : b;
}
*/

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
    
    // Rotation de Givens
    // A --> A G
    printf("A --> A G\n");
    for (int i = 0; i < n-1; i++){
        for (int j = k + i; max_int(i,j) - min_int(i,j) > 1; j--){
            printf("i = %d,j = %d, k  = %d, %d\n", i, j, k, max_int(i,j) - min_int(i,j));

            a = A[idxBand(i, j, k)];         // Élément diagonal
            b = A[idxBand(i, j+1, k)];       // Élément sous-diagonal
            printf("a = %f, b = %f\n", a, b);

            if (b == 0) {continue;}             // Si l'élément sous-diagonal est nul, on passe à l'itération suivante

            r = sqrt(a * a + b * b);
            c = a / r;
            s = -b / r;

            // Mise à jour des éléments de la matrice en appliquant la rotation
            A[idxBand(i, j, k)] = r;
            A[idxBand(i, j+1, k)] = 0;
            // On mets à jour le  reste des colonnes j et j+1
            for (int l = j + 1; l < n; l++) {//min_int(j + k + 1, n)
                A[idxBand(i, l, k)] = c * A[idxBand(i, l, k)] - s * A[idxBand(i, l + 1, k)];
                A[idxBand(i + 1, l, k)] = s * A[idxBand(i, l, k)] + c * A[idxBand(i + 1, l, k)];
                printf("i = %d, j = %d, l = %d", i, j, l);
            }
            printf("\n");
            print_band_matrix(A, n, k);


        }
        
    }
    // A G --> G^* A G = H
    for (int j = 0; j < n - 1; j++) {
        for (int i = k + j; fabs(i - j) > 2; i--) { // on commence par mettre le dernier élément de la bande à 0

            a = A[idxBand(i, j, k)];         // Élément diagonal
            b = A[idxBand(i + 1, j, k)];     // Élément sous-diagonal
        
            if (b == 0) {continue;}             // Si l'élément sous-diagonal est nul, on passe à l'itération suivante
            
            r = sqrt(a * a + b * b);
            // G^* : on prend le complexe conjugué de la matrice de Givens transposée
            c =  - a / r;
            s = b / r;
        
            // Mise à jour des éléments de la matrice en appliquant la rotation
            A[idxBand(i, j, k)] = r;
            A[idxBand(i + 1, j, k)] = 0;
            // On mets à jour le  reste des lignes i et i+1
            for (int l = j + 1; l < min_int(j + k + 1, n); l++) {
                A[idxBand(i, l, k)] = c * A[idxBand(i, l, k)] + s * A[idxBand(i, l + 1, k)];
                A[idxBand(i + 1, l, k)] = - s * A[idxBand(i, l, k)] + c * A[idxBand(i + 1, l, k)];
            }

            // Mise à jour des éléments de la matrice tridiagonale
            d[j] = A[idxBand(j, j, k)];
            e[j] = A[idxBand(j + 1, j, k)];
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

    if (m < 2) return m;
    if (square(d[m-1] + d[m-2]) + 4 * square(e[m]) < 0) return m;
    double mu = max_double((d[m-1] + d[m-2] - sqrt(square(d[m-1] + d[m-2]) + 4 * square(e[m])))/2,(d[m-1] + d[m-2] + sqrt(square(d[m] + d[m]) + 4 * square(e[m])))/2);

    for (int i = 0; i < m - 1; i++) d[i] -= mu;
    double a = d[0];
    double b = e[1];
    double c,s,r;
    for (int i = 1; i <  m; i ++){
        if (fabs(b) < eps * (fabs(a) + fabs(d[i + 1]))) continue; // L'élément est déjà nul
        r = hypot(a, b);
        c = a / r;
        s = -b / r;
        d[i - 1] = d[i - 1]*c - e[i]*s;
        e[i-1] = e[i-1]*c;
        e[i+1] = a * s + b * c;
        d[i] = c * d[i] - s * e[i];
        e[i] = s * d[i] + c * e[i];
        a = d[i];
        b = e[i + 1];
    }

    for (int i = 0; i < m - 1; i++) d[i] += mu;

    if (fabs(e[m-1]) <= eps * (fabs(d[m-2]) + fabs(d[m-1]))) {
        e[m-1] = 0.0;
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

    double *e = (double *)malloc((n-1) * sizeof(double));
    int m, index, iter = 0;

    if (e == NULL) {
        free(e);
        return -1;
    }

    tridiagonalize(A, n, k, d, e);

    for (iter = 0; iter < max_iter && n > 1 ; iter ++){
        m = step_qr_tridiag(d, e, n, eps);
        if (m  == n - 1){
            // On a isolé une valeur propre (stocké dans d), la matrice active diminue d'une taille
            // donc on incremente l'emplacement de d et e de 1
            d++;
            e++;
            index++;
            n = m;
        } else if (m != n){
            free(d);
            free(e);
            return -1;
        }

    }

    d -= index; // On ramène d à son emplacement initial -  je pense inutile ...

    free(e);
    return iter;
}