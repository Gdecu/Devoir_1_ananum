#include "devoir_1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <complex.h>


#define idxSymBand(i, j, k) ((k) * (i+1) + j ) // SYM band
#define idxSym(i, j) ((i) * ((i) + 1) / 2 + (j)) // symmetric
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

    double *A_sym = (double *)malloc(n * (n + 1) / 2 * sizeof(double));
    if (A_sym == NULL) {
        fprintf(stderr, "Erreur d'allocation mémoire pour A_copy\n");
        exit(EXIT_FAILURE);
    }

    symBand_to_sym(A, A_sym, n, k);
    //print_sym_matrix(A_sym, n, k);

    // Rotation de Givens
    // A --> G* A G = H
    for (int j = 0; j < n-1; j++){
        for (int i = n-1; i >= j + 1; i--){
            //printf("i %d, j %d\n", i, j);

            a = A_sym[idxSym(i,j)];         //A[idxSymBand(i, j, k)];         // Élément diagonal
            b = A_sym[idxSym(i+1,j)];       //A[idxSymBand(i+1, j, k)];       // Élément sous-diagonal
            //printf("a = %f, b = %f\n", a, b);


            r = sqrt(a * a + b * b);
            c = a / r;
            s = -b / r;

            if (fabs(b) < 1e-10) {
                c = 1.0; s = 0.0;
                continue;}             // Si l'élément sous-diagonal est nul, on passe à l'itération suivante

            //printf("c = %f, s = %f, r = %f\n", c, s, r);


            // Mise à jour des éléments de la matrice en appliquant la rotation

            double Aii = A_sym[idxSym(i,i)];           //A[idxSymBand(i, i, k)];
            double Ai1i = A_sym[idxSym(i+1,i)];        //A[idxSymBand(i+1, i, k)];
            double Ai1i1 = A_sym[idxSym(i+1,i+1)];     //A[idxSymBand(i+1, i+1, k)];


            double Aij, Aij1, hist;
            int col = j + 1, row = i; int reverse = 0;
            // On mets à jour le reste
            for (int x = 0; x <= n; x++){
                if (col == row) { // On ne touche pas aux éléments qui seront mis à jour plus bas
                    reverse = 1;
                    row += 2;
                    continue;
                } 
                //printf("x = %d, reverse = %d \n", x, reverse);
                if (row > n - 1) break;
                Aij = A_sym[idxSym(row,col)];                //A[idxSymBand(row, col, k)]; 
                //printf("A[%d,%d] = %f\n",row, col, Aij);
                if (reverse == 0) {
                    Aij1 = A_sym[idxSym(row+1, col)];        //A[idxSymBand(row + 1, col, k)];
                } else {
                    Aij1 = A_sym[idxSym(row,col+1)];         //A[idxSymBand(row, col + 1, k)];
                } //printf("A[%d, %d] = %f\n", row,col+1, Aij1);

                A_sym[idxSym(row,col)] = c * Aij - s * Aij1;            //A[idxSymBand(row, col, k)] = c * Aij - s * Aij1; 
                //printf("A[%d, %d], %f * %f - %f * %f = %f\n", row, col,c, Aij, s, Aij1, c * Aij - s * Aij1);

                if (reverse == 0) {
                    A_sym[idxSym(row+1,col)] = s * Aij + c * Aij1;          //A[idxSymBand(row+1, col, k)] = s * Aij + c * Aij1; 
                    //printf("A[%d, %d], %f * %f + %f * %f = %f\n", row, col + 1, s, Aij, c, Aij1, s * Aij + c * Aij1); 
                } else {
                    A_sym[idxSym(row,col+1)] = s * Aij + c * Aij1;          //A[idxSymBand(row, col + 1, k)] = s * Aij + c * Aij1; 
                    //printf("A[%d, %d], %f * %f + %f * %f = %f\n", row, col + 1, s, Aij, c, Aij1, s * Aij + c * Aij1);                   
                }


                reverse == 0 ?  col++ : row++;
            }
            reverse = 0;

            // A[i, j] = r
            A_sym[idxSym(i,j)] = (c * a - s * b);             //A[idxSymBand(i, j, k)] = c * a - s * b; 
            //printf("A[%d, %d], %f * %f - %f * %f = %f\n", i,j, c, a, s, b, A_sym[idxSym(i, j)]);

            // A[i+1, j] = 0           
            A_sym[idxSym(i+1,j)] = s * a + c * b;           //A[idxSymBand(i+1, j, k)] = s * a + c * b;
            //printf("A[%d, %d], %f * %f + %f * %f = %f\n", i+1, j, s, a, c, b, A_sym[idxSym(i+1, j)]);          
            
            // A[i, i] = c^2 * Aii - 2 * c * s * Ai1i + s^2 * Ai1i1
            A_sym[idxSym(i,i)] = c * c * Aii - 2 *c * s * Ai1i + s * s * Ai1i1;             //A[idxSymBand(i, i, k)] = c * c * Aii - 2 *c * s * Ai1i + s * s * Ai1i1;  
            //printf("A[%d, %d], %f * %f * %f - 2 * %f * %f * %f + %f * %f * %f= %f\n", i, i, c, c, Aii, c, s, Ai1i, s, s, Ai1i1, A_sym[idxSym(i, i)]);
            
            // A[i+1, i] = c * s * Aii + (c^2 - s^2) * Ai1i - c * s * Ai1i1
            A_sym[idxSym(i+1,i)] = c * s * Aii + (c*c - s*s) * Ai1i - c * s * Ai1i1;        //A[idxSymBand(i+1, i, k)] = c * s * Aii + (c*c - s*s) * Ai1i - c * s * Ai1i1; 
            
            // A[i+1, i+1] = s^2 * Aii + 2 * c * s * Ai1i + c^2 * Ai1i1
            A_sym[idxSym(i+1,i+1)] = s * s * Aii + 2 * c * s * Ai1i + c * c * Ai1i1;        //A[idxSymBand(i+1, i+1, k)] = s * s * Aii + 2 * c * s * Ai1i + c * c * Ai1i1; 


            //printf("\n");
            //print_band_matrix(A, n, k);
            //print_sym_matrix(A_sym, n, k);
        }
    }
    //print_sym_matrix(A_sym, n, k);
    for (int j = 0; j < n; j++) {
        d[j] = A_sym[idxSym(j,j)];
        if (j < n - 1) {
            e[j] = A_sym[idxSym(j + 1, j)];
        }
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
    double c,s,t_1,t_2;
    double parameters[2*m];
    t_2 = e[1];
    for (int i = 0; i <  m-1; i ++){
        c = d[i] / sqrt(square(d[i]) + square(t_2));
        s = t_2 / sqrt(square(d[i]) + square(t_2));
        parameters[2*i] = c;
        parameters[2*i+1] = s;
        d[i] = c*d[i] + s*t_2;
        t_1 = e[i+1];
        e[i+1] = s*d[i+1] + c*e[i+1];
        d[i+1] = -s*t_1 + c*d[i+1];
        if(i != m-2){
            t_2= e[i+2];
            e[i+2] *= c;
        }
    }
    for (int i = 0; i < m - 1; i++){
        c = parameters[2*i];
        s = parameters[2*i+1];
        d[i] = c * d[i] + s * e[i+1];
        e[i+1] =  s * d[i+1];
        d[i+1] *= c;    
    }
    for (int i = 0; i < m; i++) d[i] += mu;
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
int qr_eigs_(double *A, int n, int k, double eps, int max_iter, double *d) {
    double *e = (double *)malloc((n-1) * sizeof(double));
    if (e == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }
    tridiagonalize(A, n, k, d, e);
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