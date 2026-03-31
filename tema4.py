import math

def citeste_vector(nume_fisier):
    vector = []
    try:
        with open(nume_fisier, 'r') as f:
            for line in f:
                if line.strip():  
                    vector.append(float(line.strip()))
    except FileNotFoundError:
        print(f"Eroare: Fisierul {nume_fisier} nu a putut fi gasit.")
        return None
    return vector

def rezolva_sistem_rar(fisier_d0, fisier_d1, fisier_d2, fisier_b, putere_epsilon=6):
    epsilon = 10 ** (-putere_epsilon)
    
    d0 = citeste_vector(fisier_d0)
    d1 = citeste_vector(fisier_d1)
    d2 = citeste_vector(fisier_d2)
    b = citeste_vector(fisier_b)
    
    if not all([d0, d1, d2, b]):
        return 
    
    n = len(d0)
    print(f"1. Dimensiunea sistemului (n): {n}")
    
    p = n - len(d1)
    q = n - len(d2)
    print(f"2. Ordinul diagonalelor: p = {p}, q = {q}")
    
    for i, val in enumerate(d0):
        if abs(val) <= epsilon:
            print(f"3. Eroare: Element nul (sau prea mic) pe diagonala principala la indexul {i}.")
            return
    print("3. Verificare: Toate elementele diagonalei principale sunt nenule.")
    
    x = [0.0] * n  
    k = 0
    k_max = 10000
    convergenta = False
    
    while k < k_max:
        x_vechi = list(x) 
        max_delta_x = 0.0
        
        for i in range(n):
            suma = 0.0
            
            if i - p >= 0:
                suma += d1[i - p] * x[i - p]
            if i + p < n:
                suma += d1[i] * x[i + p] 
                
            if i - q >= 0:
                suma += d2[i - q] * x[i - q]
            if i + q < n:
                suma += d2[i] * x[i + q]
            
            x[i] = (b[i] - suma) / d0[i]
            
            delta_x = abs(x[i] - x_vechi[i])
            if delta_x > max_delta_x:
                max_delta_x = delta_x
                
        if max_delta_x <= epsilon:
            convergenta = True
            break
        
        if max_delta_x > 10**10:
            print(f"4. Sistemul a divergat dupa {k} iteratii (delta > 10^10).")
            return
            
        k += 1

    if convergenta:
        print(f"4. Solutia a fost aproximata cu succes dupa {k + 1} iteratii.")
        
        y = [d0[i] * x[i] for i in range(n)] 
        
        for i in range(len(d1)):
            y[i] += d1[i] * x[i + p]
            y[i + p] += d1[i] * x[i]
            
        for i in range(len(d2)):
            y[i] += d2[i] * x[i + q]
            y[i + q] += d2[i] * x[i]
            
        norma_infinit = max(abs(y[i] - b[i]) for i in range(n))
        print(f"6. Norma infinit ||A*x_GS - b||_inf : {norma_infinit}")
        

    else:
        print("Sistemul NU a convers in numarul maxim de iteratii (10.000).")

if __name__ == "__main__":
    i = 1 
    print(f"--- Rulare pentru sistemul {i} ---")
    rezolva_sistem_rar(
        fisier_d0=f"d0_{i}.txt", 
        fisier_d1=f"d1_{i}.txt", 
        fisier_d2=f"d2_{i}.txt", 
        fisier_b=f"b_{i}.txt",
        putere_epsilon=6 
    )