import numpy as np
from scipy.linalg import lu

size = 10
precision_power = 10
tolerance = 10 ** (-precision_power)

np.random.seed(1)


def generate_spd_matrix(dim):
    random_matrix = np.random.rand(dim, dim)
    return random_matrix @ random_matrix.T

def generate_responses(dim):
    return np.random.rand(dim)


def solve_with_library(matrix, responses):
    P, L, U = lu(matrix)
    solution = np.linalg.solve(matrix, responses)
    return L, U, solution

def ldlt_decomposition(mat, eps):

    n = mat.shape[0]
    diag = np.zeros(n)

    for col in range(n):

        partial_sum = 0.0
        for prev in range(col):
            partial_sum += diag[prev] * mat[col, prev] ** 2

        diag[col] = mat[col, col] - partial_sum

        if abs(diag[col]) <= eps:
            raise RuntimeError("Matricea nu e SPD")

        for row in range(col + 1, n):

            accum = 0.0
            for prev in range(col):
                accum += diag[prev] * mat[row, prev] * mat[col, prev]

            mat[row, col] = (mat[row, col] - accum) / diag[col]
    
    return diag


def solve_lower(mat, responses):

    n = len(responses)
    z = np.zeros(n)

    for i in range(n):

        s = 0.0
        for j in range(i):
            s += mat[i, j] * z[j]

        z[i] = responses[i] - s

    return z



def solve_diagonal(diagonal, vec):

    n = len(diagonal)
    y = np.zeros(n)

    for i in range(n):
        y[i] = vec[i] / diagonal[i]

    return y



def solve_upper(mat, vec):

    n = len(vec)
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):

        s = 0.0
        for j in range(i + 1, n):
            s += mat[j, i] * x[j]

        x[i] = vec[i] - s

    return x



def manual_matvec(matrix, vector):

    n = matrix.shape[0]
    result = np.zeros(n)

    for i in range(n):
        acc = 0.0
        for j in range(n):
            acc += matrix[i, j] * vector[j]
        result[i] = acc

    return result


A_original = generate_spd_matrix(size)
A_work = A_original.copy()

b = generate_responses(size)

print("Matrix A:\n", A_original)

L_lu, U_lu, x_reference = solve_with_library(A_original, b)

print("\n--- LU Decomposition ---")
print("L:\n", L_lu)
print("U:\n", U_lu)

d_vector = ldlt_decomposition(A_work, tolerance)

determinant = np.prod(d_vector)

print("\nDeterminant (LDLT):", determinant)

z_vec = solve_lower(A_work, b)
y_vec = solve_diagonal(d_vector, z_vec)
x_ldlt = solve_upper(A_work, y_vec)

print("\nSolution x_ldlt:\n", x_ldlt)

Ax = manual_matvec(A_original, x_ldlt)

Ax_b = np.linalg.norm(Ax - b)
xchol_xlib = np.linalg.norm(x_ldlt - x_reference)

print("\n||A x - b|| =", Ax_b)
print("||x_ldlt - x_lib|| =", xchol_xlib)