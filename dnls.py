# comment# Programa  DNLS

# Modulos and  import

# Add  Scrit

#  DNLS.py  +  Funcion +  nameArchivo
import sys
import random as rd
import matplotlib.pyplot as plt
import sys  # libreria  provee   variables  y funcionalidades relacionadas con el interprete
import numpy as np
from sympy import *
import cmath
import math


def RK4(h, X_00, F, nameArchivo, sigma):
    s1 = 0.9
    s2 = 4.0
    error = 0.001
    delta_t = h
    paso = 0.3
    X_N = X_00
    x_aux = X_N
    X_aux = X_N
    X_1 = X_N
    X_2 = X_N
    K_1 = F(X_N)
    k_1 = K_1
    k_2 = K_1
    k_3 = K_1
    k_4 = K_1
    K_2 = K_1
    K_3 = K_1
    K_4 = K_1
    t = 0.0

    arch = open(sigma, "w")

    arch.write(str(t))
    arch.write(" ")
    arch.write(str(0.0))
    arch.write("\n")

    arch1 = open(nameArchivo, "w")

    arch1.write(str(t))
    arch1.write(" ")

    # for i in range(len(X_N)):
    #   arch.write(str(1))
    #  arch.write(" ")

    for i in range(len(X_N)):
        arch1.write(str(X_N[i]))
        arch1.write(" ")

    arch1.write("\n")
    while t <= 20:
        K_1 = F(X_N) * delta_t
        K_2 = F(X_N + K_1 / 2.0) * delta_t
        K_3 = F(X_N + K_2 / 2.0) * delta_t
        K_4 = F(X_N + K_3) * delta_t
        X_1 = X_N + (K_1 + 2.0 * K_2 + 2.0 * K_3 + K_4) / 6.0  # siguiente paso no ms

        K_1 = F(X_N) * delta_t / 2.0
        K_2 = F(X_N + K_1 / 2.0) * delta_t / 2.0
        K_3 = F(X_N + K_2 / 2.0) * delta_t / 2.0
        K_4 = F(X_N + K_3) * delta_t / 2.0

        X_aux = X_N + (K_1 + 2.0 * K_2 + 2.0 * K_3 + K_4) * (1 / 6.0)

        K_1 = F(X_aux) * delta_t  # aqui modifique den frac de 2 a 1
        K_2 = F(X_aux + K_1 / 2.0) * delta_t
        K_3 = F(X_aux + K_2 / 2.0) * delta_t
        K_4 = F(X_aux + K_3) * delta_t
        X_2 = X_aux + (K_1 + 2.0 * K_2 + 2.0 * K_3 + K_4) / 6.0

        K_1 = F(X_N) * delta_t
        K_2 = F(X_N + 0.5 * K_1) * delta_t
        K_3 = F(X_N + 0.5 * K_2) * delta_t
        K_4 = F(X_N + K_3) * delta_t

        final = X_N + (K_1 + 2.0 * K_2 + 2.0 * K_3 + K_4) / 6.0
        X_N = final

        t = t + delta_t

        if t >= paso:
            arch.write(str(t))
            arch.write(" ")
            arch1.write(str(t))
            arch1.write(" ")

            sumanum = 0
            sumaden = 0

            for l in range(len(X_N)):
                arch1.write(str(np.absolute(X_N[l] ** 2)))
                arch1.write(" ")
                num = l * l * (np.absolute(X_N[l]) ** 2)
                den = np.absolute(X_N[l]) ** 2
                sumanum = sumanum + num
                sumaden = sumaden + den

            res = sumanum / sumaden
            arch.write(str(res))
            # arch.write(" ")

            arch.write("\n")
            arch1.write("\n")
            paso = paso + h

        else:
            continue

    #                print(t)
    return X_N


def potential(x):
    a = rd.random()
    b = rd.random()
    out = a - b
    return out


def DLS(x):
    v = 1.0 + 0.0j
    e = 0
    z = complex(0, 1)
    aux = np.zeros(101, dtype=complex)
    for i in range(len(x)):
        if i == 0:
            aux[i] = z * (v * (x[i]) + e * x[i])

        elif i == (len(x) - 1):
            aux[i] = z * (v * (x[i]) + e * x[i])

        elif i != (len(x) - 1) and i != 0:
            aux[i] = z * (v * (x[i + 1] + x[i - 1]) + e * x[i])

        # print(str(aux))
    return aux


def DLSdisorder(x):
    v = 1.0 + 0.0j
    e = 0
    z = complex(0, 1)
    aux = np.zeros(101, dtype=complex)
    for i in range(len(x)):
        if i == 0:
            aux[i] = z * (v * (x[i]) + e * x[i])

        elif i == (len(x) - 1):
            aux[i] = z * (v * (x[i]) + e * x[i])

        elif i != (len(x) - 1) and i != 0:
            aux[i] = z * (v * (x[i + 1] + x[i - 1]) + e * x[i])

        # print(str(aux))
    return aux


def DNLS(x):
    v = 1.0 + 0.0j
    e = 0
    z = complex(0, 1)
    k = [0, 1, 2, 3, 4, 5, 8, 10]
    aux = np.zeros(101, dtype=complex)
    for i in range(len(x)):
        if i == 0:
            aux[i] = z * (
                v * (x[i]) + e * x[i] + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        elif i == (len(x) - 1):
            aux[i] = z * (
                v * (x[i]) + e * x[i] + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        elif i != (len(x) - 1) and i != 0:
            aux[i] = z * (
                v * (x[i + 1] + x[i - 1])
                + e * x[i]
                + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        # print(str(aux))
    return aux


def DNLSLambda(x, params, backend=math):  #  implementar sympy
    z = complex(0, 1)
    e = params[0]
    v = params[1]
    k = params[2]
    lam = params[3]
    aux = np.zeros(101, dtype=complex)
    for i in range(len(x)):
        if i == 0:
            aux[i] = z * (
                v * (x[i]) + e * x[i] + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        elif i == (len(x) - 1):
            aux[i] = z * (
                v * (x[i]) + e * x[i] + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        elif i != (len(x) - 1) and i != 0:
            aux[i] = z * (
                v * (x[i + 1] + x[i - 1])
                + e * x[i]
                + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        # print(str(aux))
    return aux


def DNLSdisorder(x):
    v = 1.0 + 0.0j
    e = 0
    z = complex(0, 1)
    k = [0, 1, 2, 3, 4, 5, 8, 10]
    lam = 9.3
    aux = np.zeros(101, dtype=complex)
    for i in range(len(x)):
        if i == 0:
            aux[i] = z * (
                v * (x[i]) + e * x[i] + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        elif i == (len(x) - 1):
            aux[i] = z * (
                v * (x[i]) + e * x[i] + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        elif i != (len(x) - 1) and i != 0:
            aux[i] = z * (
                v * (x[i + 1] + x[i - 1])
                + e * x[i]
                + k[3] * x[i] * (np.absolute(x[i]) ** 2)
            )

        # print(str(aux))
    return aux


def mDNLS(x):
    v = 1.0 + 0.0j
    e = 0
    z = complex(0, 1)
    k = [0, 1, 2, 3, 4, 5, 8, 10]
    aux = np.zeros(101, dtype=complex)
    for i in range(len(x)):
        if i == 0:
            aux[i] = z * (
                v * (x[i])
                + e * x[i]
                + k[3]
                * x[i]
                * (2 * np.absolute(x[i]) ** 2 + np.absolute(x[i + 1]) ** 2)
            )

        elif i == (len(x) - 1):
            aux[i] = z * (
                v * (x[i])
                + e * x[i]
                + k[3] * x[i] * (2 * np.absolute(x[i]) + np.absolute(x[i - 1]) ** 2)
            )

        elif i != (len(x) - 1) and i != 0:
            aux[i] = z * (
                v * (x[i + 1] + x[i - 1])
                + e * x[i]
                + k[3]
                * x[i]
                * (
                    np.absolute(x[i + 1]) ** 2
                    + 2 * np.absolute(x[i]) ** 2
                    + np.absolute(x[i - 1]) ** 2
                )
            )

        # print(str(aux))
    return aux


def mDNLSlambda(x, params, backend=math):
    z = complex(0, 1)
    e = params[0]
    v = params[1]
    k = params[2]
    lam = params[3]
    aux = np.zeros(101, dtype=complex)
    for i in range(len(x)):
        if i == 0:
            aux[i] = z * (
                -lam * x[i]
                + v * (x[i])
                + e * x[i]
                + k[3]
                * x[i]
                * (2 * np.absolute(x[i]) ** 2 + np.absolute(x[i + 1]) ** 2)
            )
        # aux[i]=symbols('aux[i]')
        elif i == (len(x) - 1):
            aux[i] = z * (
                -lam * x[i]
                + v * (x[i])
                + e * x[i]
                + k[3] * x[i] * (2 * np.absolute(x[i]) + np.absolute(x[i - 1]) ** 2)
            )

        elif i != (len(x) - 1) and i != 0:
            aux[i] = z * (
                -lam * x[i]
                + v * (x[i + 1] + x[i - 1])
                + e * x[i]
                + k[3]
                * x[i]
                * (
                    np.absolute(x[i + 1]) ** 2
                    + 2 * np.absolute(x[i]) ** 2
                    + np.absolute(x[i - 1]) ** 2
                )
            )

        # print(str(aux))
    return aux


def condInitial(values=[1], site=[50], long=101):
    x0 = np.zeros(long)

    for i in site:
        x0[i] = values[site.index(i)]

    return x0


# parameters RK4
h = 0.01
t = 0


C0 = condInitial(values=[1, 12], site=[5, 6], long=10)

print(C0)
# p=x0.copy()
# varDNlS= RK4(h ,p, DNLSdisorder, 'DNLSdisorder.txt', 'dls24.txt');
# varmDNlS= RK4(h ,p, mDNLS, 'mDNLS.txt', 'sigmamDNLS.txt');

# var= RK4(h ,p, sys.argv[1], sys.arg[2] , 'sigma.txt');

print("\n  *** Programa Finalizado *** \n ")


#  puede ser mas adecuado  una clase de  funciones DNLS
