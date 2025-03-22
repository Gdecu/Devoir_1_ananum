# Tridiagonal matrix
import numpy as np
from scipy.linalg import solve_banded


A = np.array([
    [2, 4,  6,  0,  0, 0],
    [4, 5,  7,  9,  0, 0],
    [6, 7,  8, 10, 12, 0],
    [0, 9, 10, 11, 13, 15],
    [0, 0, 12, 13, 14, 16],
    [0, 0,  0, 15, 16, 17]])

# Givens matrix
r = np.sqrt(A[1][0] ** 2 + A[2][0] ** 2)
c = A[1][0] / r
s = -A[2][0] / r
G = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, c, -s, 0, 0, 0],
    [0, s, c, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]])

print("G = ",G)
print()
print("G.T.conj() = ",G.T.conj())
print()
print("r = ", r)
print("c = ", c)
print("s = ", s)
print(G@A)
print()
print(G@A@G.T.conj())

"""# Givens rotation for rows 2 and 3
r_23 = np.sqrt(A[2][1] ** 2 + A[3][1] ** 2)
c_23 = A[2][1] / r_23
s_23 = -A[3][1] / r_23
G_23 = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, c_23, -s_23, 0, 0],
    [0, 0, s_23, c_23, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1]])

print("G_23 = ", G_23)
print()
print("G_23.T.conj() = ", G_23.T.conj())
print()
print("r_23 = ", r_23)
print("c_23 = ", c_23)
print("s_23 = ", s_23)
print(G_23 @ A)
print()
print(G_23 @ A @ G_23.T.conj())"""