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
    
    return d, e[:-1]

# Test avec la matrice donnée
"""A = np.array([
    [  3,  6,  9, 12,  0,  0,  0],
    [  6,  7, 10, 13, 16,  0,  0],
    [  9, 10, 11, 14, 17, 20,  0],
    [ 12, 13, 14, 15, 18, 21, 25],
    [  0, 16, 17, 18, 19, 22, 24],
    [  0,  0, 20, 21, 22, 23, 26],
    [  0,  0,  0, 24, 25, 26, 27]
], dtype=float)"""

"""A = np.array([
    [2.0000, 4.0000, 6.0000, 0.0000, 0.0000, 0.0000],
    [4.0000, 5.0000, 7.0000, 9.0000, 0.0000, 0.0000],
    [6.0000, 7.0000, 8.0000, 10.0000, 12.0000, 0.0000],
    [0.0000, 9.0000, 10.0000, 11.0000, 13.0000, 15.0000],
    [0.0000, 0.0000, 12.0000, 13.0000, 14.0000, 16.0000],
    [0.0000, 0.0000, 0.0000, 15.0000, 16.0000, 17.0000]
], dtype=float)"""

# print eigvals(A)

"""A = np.array([
    [2.0000, 4.0000, 6.0000, 0.0000],
    [4.0000, 5.0000, 7.0000, 9.0000],
    [6.0000, 7.0000, 8.0000, 10.0000],
    [0.0000, 9.0000, 10.0000, 11.0000]
], dtype=float)"""

A = np.array([
    [   5.0000,   10.0000,   15.0000,   20.0000,   25.0000,   30.0000,    0.0000,    0.0000,    0.0000,    0.0000],
    [  10.0000,   11.0000,   16.0000,   21.0000,   26.0000,   31.0000,   36.0000,    0.0000,    0.0000,    0.0000],
    [  15.0000,   16.0000,   17.0000,   22.0000,   27.0000,   32.0000,   37.0000,   42.0000,    0.0000,    0.0000],
    [  20.0000,   21.0000,   22.0000,   23.0000,   28.0000,   33.0000,   38.0000,   43.0000,   48.0000,    0.0000],
    [  25.0000,   26.0000,   27.0000,   28.0000,   29.0000,   34.0000,   39.0000,   44.0000,   49.0000,   54.0000],
    [  30.0000,   31.0000,   32.0000,   33.0000,   34.0000,   35.0000,   40.0000,   45.0000,   50.0000,   55.0000],
    [   0.0000,   36.0000,   37.0000,   38.0000,   39.0000,   40.0000,   41.0000,   46.0000,   51.0000,   56.0000],
    [   0.0000,    0.0000,   42.0000,   43.0000,   44.0000,   45.0000,   46.0000,   47.0000,   52.0000,   57.0000],
    [   0.0000,    0.0000,    0.0000,   48.0000,   49.0000,   50.0000,   51.0000,   52.0000,   53.0000,   58.0000],
    [   0.0000,    0.0000,    0.0000,    0.0000,   54.0000,   55.0000,   56.0000,   57.0000,   58.0000,   59.0000]
], dtype=float)

v = np.linalg.eig(A)
d, e = tridiagonalize(A)

#print("A :\n", A)
A = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
#print("Tridiaganlize(A) :\n")
#for i in range(len(A)):
    #for j in range(len(A[0])):
        #print(f"{A[i][j]:7.2f}  ", end=" ")
    #print()
#print()
#print("eigenvalues(tridiagonalize(A)) :", np.linalg.eigvals(A))
#print("Diagonale principale :", d)
#print("Sous-diagonale :", e)

T = [
    [7.5561, 5.3277, 0.4704, 0.0000, 0.0000],
    [5.3277, 2.1896, 6.7886, 9.3469, 0.0000],
    [0.4704, 6.7886, 6.7930, 3.8350, 8.3097],
    [0.0000, 9.3469, 3.8350, 5.1942, 0.3457],
    [0.0000, 0.0000, 8.3097, 0.3457, 0.5346]
]
# Compute the eigenvalues of the tridiagonal matrix
eigenvalues = np.linalg.eigvals(T)
print("Eigenvalues of the tridiagonal matrix:", eigenvalues)
