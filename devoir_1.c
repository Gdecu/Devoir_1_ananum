#include "devoir_1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 


#define idxA(i, j) ((i)*n + (j)) // full
#define idxBand(i, j, k) ((i) * (k + 1) + (j - i + k)) // band


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
    for (int i = 0; i < n; i++) {
        d[i] = A[idxBand(i, i, k)];
        if (i < n - 1) {
            double a = A[idxBand(i, i + 1, k)];
            double b = A[idxBand(i + 1, i, k)];
            double r = hypot(a, b);
            double c = a / r;
            double s = -b / r;
            
            A[idxBand(i, i + 1, k)] = r;
            e[i] = r;
            
            for (int j = i + 1; j < n; j++) {
                if (j - i <= k) {
                    double A_ji = A[idxBand(j, i, k)];
                    double A_j1i = (j + 1 < n) ? A[idxBand(j + 1, i, k)] : 0.0;
                    
                    A[idxBand(j, i, k)] = c * A_ji - s * A_j1i;
                    if (j + 1 < n) {
                        A[idxBand(j + 1, i, k)] = s * A_ji + c * A_j1i;
                    }
                }
            }
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
    return 0;
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