#include "devoir_1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 

#define idxA(i, j) ((i)*n + (j)) // full
#define idxBand(i, j, k) ((i) * (k + 1) + (j - i + k)) // band


/**
 * @brief calcule le carré d'un nombre
 * @param x nombre dont on veut le carré
 * @return le carré de x
 */
double square( double x){
    return x * x;
}



/**
 * @brief calcule le maximum entre deux entiers
 * @param a premier entier
 * @param b deuxième entier
*/
double max_double(double a, double b) {
    return (a > b) ? a : b;
}

// Applique une rotation de Givens pour annuler A[i+1, i]
void givens_rotation(double *A, int n, int k, int i, double *c, double *s) {
    double a = A[idxBand(i, i, k)];   // Élément diagonal
    double b = A[idxBand(i + 1, i, k)]; // Élément sous-diagonal
    
    double r = sqrt(a * a + b * b);
    *c = a / r;
    *s = -b / r;
    
    // Mise à jour des éléments de la matrice en appliquant la rotation
    A[idxBand(i, i, k)] = r;
    A[idxBand(i + 1, i, k)] = 0;  // Annulation de l'élément sous-diagonal
}

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
    double c, s;
    
    for (int i = 0; i < n - 1; i++) {
        givens_rotation(A, n, k, i, &c, &s);
        
        // Stockage des valeurs tridiagonales
        d[i] = A[idxBand(i, i, k)];
        e[i] = (i < n - 1) ? A[idxBand(i + 1, i, k)] : 0;
    }
    d[n - 1] = A[idxBand(n - 1, n - 1, k)]; // Dernière diagonale
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

    if (m == 1){
        return 1;
    }

    if (square(d[m-1] + d[m-2]) + 4 * square(e[m-1]) < 0){
        return m;
    }

    double mu = max_double((d[m-1] + d[m-2] - sqrt(square(d[m-1] + d[m-2]) + 4 * square(e[m-1])))/2,(d[m-1] + d[m-2] + sqrt(square(d[m-1] + d[m-2]) + 4 * square(e[m-1])))/2);

    for (int i = 0; i < m - 1; i++){
        d[i] -= mu;
    }

    for (int i = 0; i < m - 1; i++){
        d[i] += mu;
    }

    
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