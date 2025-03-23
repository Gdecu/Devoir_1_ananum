# Tridiagonal matrix
import numpy as np
from scipy.linalg import solve_banded


A = np.array([
    [  3,  6,  9, 12,  0,  0,  0],
    [  6,  7, 10, 13, 16,  0,  0],
    [  9, 10, 11, 14, 17, 20,  0],
    [ 12, 13, 14, 15, 18, 21, 25],
    [  0, 16, 17, 18, 19, 22, 24],
    [  0,  0, 20, 21, 22, 23, 26],
    [  0,  0,  0, 24, 25, 26, 27]
], dtype=float)


# Givens matrix
r = np.sqrt(A[2][0] ** 2 + A[3][0] ** 2)
c = A[2][0] / r
s = -A[3][0] / r
G = np.array([
    [1, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [0, 0, c, -s, 0, 0, 0],
    [0, 0, s, c, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 1]])

print("A = \n", A)
print()
print("r = ", r," c = ", c," s = ", s)
print()
H = G@A@G.T.conj()
for i in range(7):
    for j in range(7):
        print(f"{H[i][j]:7.2f}", end=" ")
    print()
print()

# Givens rotation for rows 2 and 3
r_23 = np.sqrt(A[2][1] ** 2 + A[3][1] ** 2)
c_23 = A[2][1] / r_23
s_23 = -A[3][1] / r_23
G_23 = np.array([
    [1, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [0, 0, c_23, -s_23, 0, 0, 0],
    [0, 0, s_23, c_23, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 1]])


print("r_23 = ", r_23, "c_23 = ", c_23, "s_23 = ", s_23)
H = G_23@A@G_23.T.conj()
for i in range(7):
    for j in range(7):
        print(f"{H[i][j]:7.2f}", end=" ")
    print()
print()


# Givens rotation for rows 2 and 4
r_24 = np.sqrt(A[2][1] ** 2 + A[4][1] ** 2)
c_24 = A[2][1] / r_24
s_24 = -A[4][1] / r_24

G_24 = np.array([
    [1, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [0, 0, c_24, -s_24, 0, 0, 0],
    [0, 0, s_24, c_24, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 1]])

print("r_24 = ", r_24, "c_24 = ", c_24, "s_24 = ", s_24)
H = G_24@A@G_24.T.conj()
for i in range(7):
    for j in range(7):
        print(f"{H[i][j]:7.2f}", end=" ")
    print()
print()