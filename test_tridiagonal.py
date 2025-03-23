import numpy as np

def givens_rotation(a, b):
    """Compute the Givens rotation matrix elements c and s such that:
    G @ [a, b]ᵀ = [r, 0]ᵀ
    """
    if b == 0:
        return 1, 0
    elif abs(b) > abs(a):
        tau = -a / b
        s = 1 / np.sqrt(1 + tau**2)
        c = s * tau
    else:
        tau = -b / a
        c = 1 / np.sqrt(1 + tau**2)
        s = c * tau
    return c, s

def tridiagonalize(A):
    """Transform a symmetric band matrix A into a tridiagonal matrix using Givens rotations."""
    n = A.shape[0]
    A = A.copy()
    d = np.zeros(n)
    e = np.zeros(n)
    
    for j in range(n - 2):
        for i in range(j + 2, n):
            if A[i, j] != 0:
                c, s = givens_rotation(A[j + 1, j], A[i, j])
                
                # Apply rotation on rows (left multiplication)
                for k in range(j, n):
                    A_jk, A_ik = A[j + 1, k], A[i, k]
                    A[j + 1, k] = c * A_jk - s * A_ik
                    A[i, k] = s * A_jk + c * A_ik
                
                # Apply rotation on columns (right multiplication)
                for k in range(n):
                    A_kj, A_ki = A[k, j + 1], A[k, i]
                    A[k, j + 1] = c * A_kj - s * A_ki
                    A[k, i] = s * A_kj + c * A_ki
    
    # Extract diagonal and sub-diagonal
    d[:] = np.diag(A)
    e[:-1] = np.diag(A, k=1)
    
    return d, e

# Test avec la matrice donnée
A = np.array([
    [4.5865   , 0.4704 ,   9.3469  ,  8.3097  ,  0.0000 ],
   [0.4704  ,  6.7886   , 3.8350  ,  0.3457   , 6.7115 ],
   [9.3469   , 3.8350  ,  5.1942  ,  0.5346  ,  0.0770 ],
   [8.3097  ,  0.3457  ,  0.5346  ,  5.2970   , 3.8342 ],
   [0.0000  ,  6.7115   , 0.0770   , 3.8342   , 0.6684]
], dtype=float)

v = np.linalg.eig(A)
d, e = tridiagonalize(A)

print("Valeur propre de A :", v[0])
#print("Diagonale principale :", d)
#print("Sous-diagonale :", e)
