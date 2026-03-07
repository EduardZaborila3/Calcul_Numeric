import math
import random
import time


def precizia_masina():
    u = 1.0

    while 1.0 + u != 1.0:
        u /= 10.0

    return u * 10


u = precizia_masina()

print("Precizia masina u =", u)



x = 1.0
y = u / 10
z = u / 10

print("(x + y) + z =", (x + y) + z)
print("x + (y + z) =", x + (y + z))
print("Sunt egale?", ((x + y) + z) == (x + (y + z)))



x = 1e200
y = 1e-200
z = 15

print("(x * y) * z =", (x * y) * z)
print("x * (y * z) =", x * (y * z))
print("Sunt egale?", ((x * y) * z) == (x * (y * z)))




def my_tan_lentz(x, eps=1e-12):

    mic = 1e-12

    x = ((x + math.pi/2) % math.pi) - math.pi/2

    if abs(abs(x) - math.pi/2) < 1e-10:
        return float('inf')

    b0 = 1.0
    f = b0 if b0 != 0 else mic

    C = f
    D = 0.0
    j = 1

    while True:
        a = -x * x
        b = 2 * j + 1

        D = b + a * D
        if D == 0:
            D = mic
        D = 1.0 / D

        C = b + a / C
        if C == 0:
            C = mic

        delta = C * D
        f *= delta

        if abs(delta - 1.0) < eps:
            break

        j += 1

    return x / f


def my_tan_poly(x):

    x = ((x + math.pi/2) % math.pi) - math.pi/2

    if abs(abs(x) - math.pi/2) < 1e-10:
        return float('inf')

    if abs(x) > math.pi/4:
        if x > 0:
            return 1 / my_tan_poly(math.pi/2 - x)
        else:
            return -1 / my_tan_poly(math.pi/2 + x)

    c1 = 0.33333333333333333
    c2 = 0.133333333333333333
    c3 = 0.053968253968254
    c4 = 0.0218694885361552

    x2 = x * x
    x3 = x2 * x
    x4 = x2 * x2
    x6 = x4 * x2

    return x + x3*(c1 + x2 * c2 + x4 * c3 + x6 * c4) # x + c1 * x3 + c2 * x5 + c3 * x7 + c4 * x9




N = 10000
values = [random.uniform(-math.pi/2, math.pi/2) for _ in range(N)]

#LENTZ
start = time.time()
errors_lentz = []

for x in values:
    real = math.tan(x)
    approx = my_tan_lentz(x)
    errors_lentz.append(abs(real - approx))

time_lentz = time.time() - start

print("LENTZ:")
print("Eroare medie:", sum(errors_lentz)/N)
print("Timp executie:", time_lentz, "secunde\n")

#POLINOM
start = time.time()
errors_poly = []

for x in values:
    real = math.tan(x)
    approx = my_tan_poly(x)
    errors_poly.append(abs(real - approx))

time_poly = time.time() - start


print("POLINOM:")
print("Eroare medie:", sum(errors_poly)/N)
print("Timp executie:", time_poly, "secunde\n")