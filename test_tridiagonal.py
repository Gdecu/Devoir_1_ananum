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
A = np.array([
    [  3,  6,  9, 12,  0,  0,  0],
    [  6,  7, 10, 13, 16,  0,  0],
    [  9, 10, 11, 14, 17, 20,  0],
    [ 12, 13, 14, 15, 18, 21, 25],
    [  0, 16, 17, 18, 19, 22, 24],
    [  0,  0, 20, 21, 22, 23, 26],
    [  0,  0,  0, 24, 25, 26, 27]
], dtype=float)

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

d, e = tridiagonalize(A)
# print eigenvalues(tridiagonalize(A))
print("A :\n", A)
A = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
print("Tridiaganlize(A) :\n")
for i in range(len(A)):
    for j in range(len(A[0])):
        print(f"{A[i][j]:7.2f}  ", end=" ")
    print()
print()
#print("eigenvalues(tridiagonalize(A)) :", np.linalg.eigvals(A))
print("Diagonale principale :", d)
print("Sous-diagonale :", e)

# Construct the tridiagonal matrix from the diagonal and sub-diagonal
#T = np.diag([ 2.0000,  13.5385,  24.0247,  18.3403,  -1.6067,  18.3403], k = 0) + np.diag([7.2111,  16.6923,  22.4375,   0.5719,  24.0247], k = -1) + np.diag([7.2111,  16.6923,  22.4375,   0.5719,  24.0247], k = 1)
#print("Tridiagonal matrix:\n", T)
# Compute the eigenvalues of the tridiagonal matrix
#eigenvalues = np.linalg.eigvals(T)
#print("Eigenvalues of the tridiagonal matrix:", eigenvalues)
