import numpy as np

def householder_qr(A_init, epsilon=1e-12):
    n = A_init.shape[0]
    A = A_init.copy()
    Q_bar = np.eye(n)
    
    for r in range(n - 1):
        sigma = np.sum(A[r:, r]**2)
        
        #verific daca matricea e singulara
        if sigma <= epsilon:
            continue
            
        k = np.sqrt(sigma)
        if A[r, r] > 0:
            k = -k
            
        beta = sigma - k * A[r, r]
        u = np.zeros(n)
        u[r] = A[r, r] - k
        u[r+1:] = A[r+1:, r]
        
        #parcurg toate coloanele ramase la dreapta
        for j in range(r + 1, n):
            suma = 0.0
            for i in range(r, n):
                suma += u[i] * A[i, j]
            gamma = suma / beta
            for i in range(r, n):
                A[i, j] = A[i, j] - gamma * u[i]
            
        A[r, r] = k
        A[r+1:, r] = 0
        
        for j in range(n):
            suma_q = 0.0
            for i in range(r, n):
                suma_q += u[i] * Q_bar[i, j]
            gamma = suma_q / beta
            for i in range(r, n):
                Q_bar[i, j] = Q_bar[i, j] - gamma * u[i]
            
    R = A
    Q = Q_bar.T
    return Q, R

#rezolvare Rx=b
def solve_upper_triangular(R, b):
    n = len(b)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(R[i, i+1:], x[i+1:])) / R[i, i]
    return x

def get_inverse(Q, R):
    n = Q.shape[0]
    A_inv = np.zeros((n, n))
    for j in range(n):
        e_j = np.zeros(n)
        e_j[j] = 1
        b = Q.T @ e_j 
        #rezolv Rx=b
        x_star = solve_upper_triangular(R, b)
        A_inv[:, j] = x_star
    return A_inv

def main():
    n = 100
    print(f"Rulare pentru n = {n}")
    A_init = np.random.rand(n, n) * 10
    s = np.random.rand(n) * 10
    
    #calculez b
    b_init = A_init @ s 
    
    #descompunereaQR
    Q_h, R_h = householder_qr(A_init)
    #rezolvare sistem cu alg meu
    b_h = Q_h.T @ b_init
    x_householder = solve_upper_triangular(R_h, b_h)
    
    #rezolvare sistem cu biblioteca
    Q_lib, R_lib = np.linalg.qr(A_init)
    b_lib = Q_lib.T @ b_init
    x_qr = solve_upper_triangular(R_lib, b_lib)
    
    norma_solutii = np.linalg.norm(x_qr - x_householder, 2)
    print(f"x_QR - x_Householder = {norma_solutii:.10e}")
    
    #erori
    err1 = np.linalg.norm(A_init @ x_householder - b_init, 2)
    err2 = np.linalg.norm(A_init @ x_qr - b_init, 2)
    
    norm_s = np.linalg.norm(s, 2)
    err3 = np.linalg.norm(x_householder - s, 2) / norm_s
    err4 = np.linalg.norm(x_qr - s, 2) / norm_s
    
    print("\nErori:")
    print(f" - A*x_Householder - b = {err1:.10e}")
    print(f" - A*x_QR - b = {err2:.10e}")
    print(f" - x_Householder - s / s = {err3:.10e}")
    print(f" - x_QR - s / s= {err4:.10e}")
    
    #inversa matr A
    A_inv_householder = get_inverse(Q_h, R_h)
    A_inv_lib = np.linalg.inv(A_init)
    
    norma_inverselor = np.linalg.norm(A_inv_householder - A_inv_lib, 2)
    print("\n Inversa matricei A:")
    print(f" - A_inv_Householder - A_inv_bibl = {norma_inverselor:.10e}")

if __name__ == "__main__":
    main()