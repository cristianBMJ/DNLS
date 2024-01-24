from pyneqsys.symbolic import *
import numpy as np


import matplotlib.pyplot as plt
import math
from numpy import linalg as LA


# function Discrete Non-linear Schrodinger
def DNLS(x, P, backend=math):
    # include information help
    x = np.array(x)
    P = np.array(P)
    N = len(x) - 1
    func = lambda x, P: [
        (P[0] - P[3]) * x[n]
        + P[1] * (x[n - 1] + x[n + 1])
        + P[2]
        * x[n]
        * (
            2 * np.absolute(x[n]) ** 2
            + np.absolute(x[n - 1]) ** 2
            + np.absolute(x[n + 1]) ** 2
        )
        for n in range(len(x))
        if (n != 0 and n != len(x) - 1)
    ]
    result = np.array(func(x, P))
    result = np.append(
        result,
        (P[0] - P[3]) * x[N]
        + P[1] * (x[N - 1])
        + P[2] * x[N] * (2 * np.absolute(x[N]) ** 2 + np.absolute(x[N - 1]) ** 2),
    )
    result = np.insert(
        result,
        0,
        (P[0] - P[3]) * x[0]
        + P[1] * (x[1])
        + P[2] * x[0] * (2 * np.absolute(x[0]) ** 2 + np.absolute(x[1]) ** 2),
    )
    return result


# =============================================================================
# condiciones iniciales
xx0 = np.zeros(101)
xx0[0] = 1
# x0[1]=1

y0 = np.zeros(101)
y0[50] = 0
y0[51] = 1
w0 = np.zeros(101)
w0[50] = 1
w0[51] = -1
z0 = np.zeros(101)
z0[49] = 1
z0[50] = 1
z0[51] = -1
z0[52] = -1

P0 = [0, 1, 2, 9.3]  # [ eSite, V , X , lam ]
print("potencia es: ", P0)
N = len(xx0)
Np = len(P0)
DNLS_sys = SymbolicSys.from_callback(DNLS, N, Np)

# Method Newton, to solved equation non-linear
xx, infox = DNLS_sys.solve(xx0, P0)
site = np.arange(0, len(xx) - 1, 1)
y, infoy = DNLS_sys.solve(y0, P0)
w, infow = DNLS_sys.solve(w0, P0)
z, infoz = DNLS_sys.solve(z0, P0)
# Grafica
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(site[45:55], xx[45:55])
axs[0, 1].plot(site[45:55], y[45:55])
axs[1, 0].plot(site[45:55], w[45:55])
axs[1, 1].plot(site[40:60], z[40:60])
plt.show()

# =============================================================================
# Modes Surface
# =============================================================================
# =============================================================================
# x0=np.zeros(101)
# x0[0]=1
# x0[1]=0
# y0=np.zeros(101)
# y0[0]=1
# y0[1]=1
# w0=np.zeros(101)
# w0[0]=1
# w0[1]=-1
# z0=np.zeros(101)
# z0[0]=1
# z0[1]=1
# z0[2]=-1
# z0[3]=-1
#
# N=len(x0)
# Np=len(P)
# DNLS_sys = SymbolicSys.from_callback(DNLS, N, Np)
# x, info = DNLS_sys.solve(x0,P)
# site=np.arange(0,len(x)-1,1)
# y, info = DNLS_sys.solve(y0,P)
# w, info = DNLS_sys.solve(w0,P)
# z, info = DNLS_sys.solve(z0,P)
# #Grafica
# fig,axs=plt.subplots(2,2)
# axs[0,0].plot( site[0:10] ,x[0:10])
# axs[0,1].plot( site[0:10] ,y[0:10])
# axs[1,0].plot( site[0:10] ,w[0:10])
# axs[1,1].plot( site[0:15] ,z[0:15])
# plt.xlabel('Site')
# plt.show()
#
# # =============================================================================
# =============================================================================
# Matrices A and B
# =============================================================================

kron = lambda i, j: 1 if (i == j) else 0
conj = lambda z: np.conjugate(z)


def matrices(q, P):
    V = P[1]
    chi = P[2]
    lam = P[3]
    A = np.zeros((len(q), len(q)))
    B = np.zeros((len(q), len(q)))
    N = len(q) - 1
    for n in range(1, N):
        for m in range(len(q)):
            A[n][m] = (
                (
                    -lam
                    + chi
                    * (
                        np.absolute(q[n + 1]) ** 2
                        + np.absolute(q[n - 1]) ** 2
                        + 4 * np.absolute(q[n]) ** 2
                        - 2 * chi * q[n] ** 2
                    )
                )
                * kron(n, m)
                + (V + chi * q[n] * (conj(q[n + 1]) - q[n + 1])) * kron(n, m - 1)
                + (V + chi * q[n] * (conj(q[n - 1]) - q[n - 1])) * kron(n, m + 1)
            )
            B[n][m] = (
                (
                    lam
                    - chi
                    * (
                        np.absolute(q[n + 1]) ** 2
                        + np.absolute(q[n - 1]) ** 2
                        + 4 * np.absolute(q[n]) ** 2
                        - 2 * chi * q[n] ** 2
                    )
                )
                * kron(n, m)
                - (V + chi * q[n] * (q[n + 1] + q[n + 1])) * kron(n, m - 1)
                - (V + chi * q[n] * (q[n - 1] + q[n - 1])) * kron(n, m + 1)
            )

    A[0][1] = V + chi * q[0] * (conj(q[1]) - q[1])
    B[0][1] = -(V + chi * q[0] * (conj(q[1]) + q[1]))
    A[N][N - 1] = V + chi * q[N] * (conj(q[N - 1]) - q[N - 1])
    B[N][N - 1] = -(V + chi * q[N] * (conj(q[N - 1]) + q[N - 1]))

    A[0][0] = -lam + chi * (
        np.absolute(q[1]) ** 2 + 4 * np.absolute(q[0]) ** 2 - 2 * chi * q[0] ** 2
    )
    A[N][N] = -lam + chi * (
        np.absolute(q[N - 1]) ** 2 + 4 * np.absolute(q[N]) ** 2 - 2 * chi * q[N] ** 2
    )
    B[0][0] = lam - chi * (
        np.absolute(q[1]) ** 2 + 4 * np.absolute(q[0]) ** 2 - 2 * chi * q[0] ** 2
    )
    B[N][N] = lam - chi * (
        np.absolute(q[N - 1]) ** 2 + 4 * np.absolute(q[N]) ** 2 - 2 * chi * q[N] ** 2
    )

    W, V = LA.eig(B * A)
    return np.sort((np.absolute(W)))  # probar logaritmo


def potencia(g, X0, p):
    N = len(X0)
    Np = len(p)
    pot = np.array([])
    lam = np.array([])
    alpha = 1
    h = 0.1
    DNLS_sys = SymbolicSys.from_callback(DNLS, N, Np)
    phi0, info = DNLS_sys.solve(X0, p)
    for L in g:
        # phi=0
        aux = 0
        # =============================================================================
        #         u=np.sqrt(L)
        #         eps=np.sinh(u*h)/u
        #         nn=-2*np.cosh(u*h) + alpha*eps
        #         p[3]=-nn+2
        # =============================================================================
        p[3] = np.sqrt(L)
        # p[3]= L
        print(P[3], "\n")
        phi, info = DNLS_sys.solve(phi0, p)
        for i in phi:
            aux += np.absolute(i) ** 2
        # =============================================================================
        #         if  aux!=0:
        #             pot=np.append(pot,aux)
        #             lam=np.append(lam,L)
        # =============================================================================

        pot = np.append(pot, aux)
        lam = np.append(lam, L / 10)
        phi0 = phi

    # p[3]=9.3
    return lam, pot


def potenciaSurface(g):
    res = np.array([])
    for lam in g:
        # A=lam -1/lam # aproximacion
        A = lam / 2 - np.sqrt(lam**2 / 4 - 1)
        aux = 0
        for n in range(len(g)):
            aux += A * np.e ** (-4 * n * np.log(np.sqrt(A)))
        res = np.append(res, aux)

    return g, res


def potenciaBulk(g):  # encoding
    res = np.array([])
    for lam in g:
        A = lam - 1 / lam
        aux = 0
        for n in range(len(g)):
            aux += A * np.e ** (-4 * n * np.log(np.sqrt(A)))
        res = np.append(res, aux)

    return g, res


# ===================
# =============================================================================
#
# def autopotencia( x0 , P ):
#     N=len(x0)
#     Np=len(P)
#     pot=np.array([])
#     values=np.array([])
#     Lam=P[3]
#     x_aux=x0
#     for i in range(N):
#         aux=0
#         P[3]=Lam
#         g=matrices( x_aux , 1, 1 , Lam )
#         values=np.append(values, g[0])
#         DNLS_sys = SymbolicSys.from_callback(DNLS, N, Np)
#         x, info = DNLS_sys.solve(x_aux,P)
#         for i in x:
#             aux+=np.absolute(i)**2
#         pot=np.append(pot,aux)
#         Lam=g[0]
#         x_aux=x
#
#     return values,pot
# =============================================================================

P = P0.copy()
# surface

Y = y.copy()
g0 = matrices(Y, P)
# g=matrices(x,P)
g = np.linspace(2.1, 8, 100)
# =============================================================================
#
Lam, Power = potencia(g0, y0, P)
# LamS, PowerS= potenciaSurface(g0 )
# #Lam, Power= potenciaSurface(g)
fig2 = plt.plot(Lam, Power, "*")
# #fig3=plt.plot(LamS, PowerS)
plt.show()
#
# =============================================================================
print(P0)
print("\n ****Finalizado!**** \n")
## tiene -2 y  lambda = 5.5,  bonito
