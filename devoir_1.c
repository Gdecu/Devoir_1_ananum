#include "devoir_1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <complex.h>


#define idxSymBand(i, j, k) ((k) * (i+1) + j ) 
#define idxSym(i, j) ((i) * ((i) + 1) / 2 + (j)) 
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
void tridiagonalize_band(double *A, int n, int k, double *d, double *e) {
     double c, s, a, b, r;
     double *A_sym = (double *)malloc(n * (n + 1) / 2 * sizeof(double));
     if (A_sym == NULL) {
         fprintf(stderr, "Erreur d'allocation mémoire pour A_copy\n");
         exit(EXIT_FAILURE);
     }
     symBand_to_sym(A, A_sym, n, k);
     for (int j = 0; j < n-1; j++){
         for (int i = n-1; i >= j + 1; i--){
             a = A_sym[idxSym(i,j)];    
             b = A_sym[idxSym(i+1,j)];      
             r = sqrt(a * a + b * b);
             c = a / r;
             s = -b / r;
             if (fabs(b) < 1e-10) {
                 c = 1.0; s = 0.0;
                 continue;}         
             double Aii = A_sym[idxSym(i,i)];         
             double Ai1i = A_sym[idxSym(i+1,i)];    
             double Ai1i1 = A_sym[idxSym(i+1,i+1)];    
             double Aij, Aij1, hist;
             int col = j + 1, row = i; int reverse = 0;
             for (int x = 0; x <= n; x++){
                 if (col == row) { 
                     reverse = 1;
                     row += 2;
                     continue;
                 } 
                 if (row > n - 1) break;
                 Aij = A_sym[idxSym(row,col)];             
                 if (reverse == 0) {
                     Aij1 = A_sym[idxSym(row+1, col)];      
                 } else {
                     Aij1 = A_sym[idxSym(row,col+1)];        
                 } 
                 A_sym[idxSym(row,col)] = c * Aij - s * Aij1;            
                 if (reverse == 0) {
                     A_sym[idxSym(row+1,col)] = s * Aij + c * Aij1;          
                 } else {
                     A_sym[idxSym(row,col+1)] = s * Aij + c * Aij1;                           
                 }
                 reverse == 0 ?  col++ : row++;
             }
             reverse = 0;
             A_sym[idxSym(i,j)] = (c * a - s * b);    
             A_sym[idxSym(i+1,j)] = s * a + c * b;       
             A_sym[idxSym(i,i)] = c * c * Aii - 2 *c * s * Ai1i + s * s * Ai1i1;  
             A_sym[idxSym(i+1,i)] = c * s * Aii + (c*c - s*s) * Ai1i - c * s * Ai1i1;   
             A_sym[idxSym(i+1,i+1)] = s * s * Aii + 2 * c * s * Ai1i + c * c * Ai1i1;       
         }
     }
     for (int j = 0; j < n; j++) {
         d[j] = A_sym[idxSym(j,j)];
         if (j < n - 1) {
             e[j] = A_sym[idxSym(j + 1, j)];
         }
     }
     for (int i = n - 1; i > 0; i--) {
            e[i] = e[i - 1];
     }
     free(A_sym);
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
     double mu;  
     double delta = (d[m - 2]-d[m - 1]) / 2.0;
     double sign_delta = (delta >= 0) ? 1.0 : -1.0;
     if(delta <= 10e-12){
         mu = d[m - 1] - fabs(e[m - 1]);
     }else{
         mu = d[m - 1] - ((sign_delta*square(e[m - 1]))) / (fabs(delta) +  sqrt(square(delta) + square(e[m - 1])));
     }
     for (int i = 0; i < m; i++) d[i] -= mu;
     double c,s,r,t_1,t_2;
     double cosinus[m];
     double sinus[m];    
     t_2 = e[1];
     for (int i = 0; i <  m-1; i ++){
         r = sqrt(square(d[i]) + square(t_2));
         c = d[i] / r;
         s = t_2 / r;
         cosinus[i] = c;
         sinus[i] = s;
         d[i] = c * d[i] + s * t_2;
         t_1 = e[i+1];
         e[i+1] = s * d[i+1] + c * e[i+1];
         d[i+1] = -s * t_1 + c * d[i+1];
         if(i < m-2){
             t_2= e[i+2];
             e[i+2] *= c;
         }
     }
     for (int i = 0; i < m - 1; i++){
         c = cosinus[i];
         s = sinus[i];
         d[i] = c * d[i] + s * e[i+1];
         e[i+1] =  s * d[i+1];
         d[i+1] *= c;    
     }
     for (int i = 0; i < m; i++) d[i] += mu;
     e[0] = 0.0;
     for(int i = m-1; i>=0; i--){
         if( fabs(e[i]) < eps * (fabs(d[i-1]) + fabs(d[i])) ) {
             m--;
         }
         else{break;}
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
int qr_eigs_band(double *A, int n, int k, double eps, int max_iter, double *d) {
    double *e = (double *)malloc((n-1) * sizeof(double));
    if (e == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }
    tridiagonalize_band(A, n, k, d, e);
    for (int i = 0; i < max_iter; i++) {
         int m = step_qr_tridiag(d, e, n, eps);
         if (m == 0) {
            return i;
         }
         n = m;
    }
    free(e);
    return -1;
}