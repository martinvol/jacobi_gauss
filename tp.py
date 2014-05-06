# -*- coding: utf-8 -*-

# Integrantes:
# Gonzalo Guzzardi (94258)
# Martín Volpe     (95442)

# Ejecutar bajo Python 2.x

import csv
import numpy as np
from numpy.core import Inf
from numpy.linalg import norm
from math import log, exp

B = np.matrix(
            [ [2    ,   -1./2,   0     ],
              [-1./2,   2    ,   -1./2 ],
              [0    ,   -1./2,   2     ],
            ],
            dtype=np.float64,
)
C = np.matrix(
            [ [4 ,  -1,   0  ],
              [-1,  4 ,   -1 ],
              [0 ,  -1,   4  ],
            ],
            dtype=np.float64,
)
D = np.matrix(
            [ [3    ,  -1./2,   0     ],
              [-1./2,  3    ,   -1./2 ],
              [0    ,  -1./2,   3     ],
            ],
            dtype=np.float64,
)

E = np.matrix(
            [ [-1,  0 ,  0  ],
              [0 ,  -1,  0  ],
              [0 ,  0 ,  -1 ],
            ],
            dtype=np.float64,
)

A6 = np.concatenate((np.concatenate((B,E), axis=1),
                np.concatenate((E,D), axis=1)), axis=0)


VECTOR_ORIGEN = np.array([0,0,0], dtype=np.float64,)

def crear_matriz(n):
    """Crea una matriz con las especificaciones del enunciado"""
    matriz = np.concatenate((B,E), axis=1)
    matriz = np.concatenate((matriz, np.zeros((3,n-6))), axis=1)

    if n > 6:
        cosnt_medio = np.concatenate((E,C,E), axis=1)
        for i in xrange(0, n-6, 3):
            medio = np.concatenate((np.zeros((3,i)), cosnt_medio), axis=1) if i > 0 else cosnt_medio
            medio = np.concatenate((medio, np.zeros((3,n -(9+i)))), axis=1) if n -(9+i) > 0 else medio
            matriz = np.concatenate((matriz, medio), axis=0)
    
    fin  = np.concatenate((E,D), axis=1)
    fin = np.concatenate((np.zeros((3,n-6)), fin), axis=1)

    matriz = np.concatenate((matriz, fin), axis=0)

    return matriz


def jacobi_paso(matrix, b, x=None):
    """Realiza una iteración con el método de Jacobi"""

    x_temp = np.array(x, dtype=np.float64)

    for i, fila in enumerate(matrix):
        restar = 0
        for j in range(1, fila.size):
            indice = (j+i)%fila.size # estos son los subindices a encontrar

            if 0 <= i < 3 and indice < 6:
                # no incluimos los ceros del principio de la matriz
                restar -= x.item(indice)*fila.item(indice)
            elif fila.size -1 -3 < i < fila.size -1 and indice > fila.size -1 -6:
                # no incluimos los ceros del final de la matriz
                restar -= x.item(indice)*fila.item(indice)
            elif not (i%3 == 0 and j%3 == 2) or (i%3 == 2 and j%3 == 0):
                restar -= x.item(indice)*fila.item(indice)

        temp = b.item(i) + restar
        x_temp[i] = temp/float(fila.item(i))
        # print x_temp # descomentar esta línea para ver las iteraciones
    return x_temp

def gs_paso(matrix, b, x=None):
    """Realiza una iteración con el método de Gauss-Seidel"""

    x_temp = np.array(x, dtype=np.float64)

    for i, fila in enumerate(matrix):
        restar = 0
        for j in range(1, fila.size):
            indice = (j+i)%fila.size # estos son los subindices a encontrar

            if 0 <= i < 3 and indice < 6:
                # no incluimos los ceros del principio de la matriz
                restar -= x_temp.item(indice)*fila.item(indice)
            elif fila.size -1 -3 < i < fila.size -1 and indice > fila.size -1 -6:
                # no incluimos los ceros del final de la matriz
                restar -= x_temp.item(indice)*fila.item(indice)
            elif not (i%3 == 0 and j%3 == 2) or (i%3 == 2 and j%3 == 0):
                restar -= x_temp.item(indice)*fila.item(indice)

        temp = b.item(i) + restar
        x_temp[i] = temp/float(fila.item(i))
        # print x_temp # descomentar esta línea para ver las iteraciones

    return x_temp


def rtol(xk, xk_1):
    """Calcula el error entre los vectores de dos iteraciones"""
    return norm(xk - xk_1 , ord=Inf)/norm(xk, ord=Inf)

def radio_espectral(xk, x, x0, k):
    """Calcula el radio espectral"""
    a = norm(xk - x, ord=Inf)
    b = norm(x0 - x, ord=Inf)
    lnp = (log(a) - log(b))/k
    return exp(lnp)

def resolver(algoritmo, matriz, b, x):
    """resuelve un sistema de ecuaciones y retorna sus radios espectrales
    y errores para cada iteración."""
    radios_espectrales = {}
    errores = {}
    # k = 0, # numero de iteraciones
    x0 = np.zeros(matriz[0].size)
    xk_1 = algoritmo(matriz, b, x0)
    # k = 1
    
    errores[1] = rtol(xk_1, x0)
    radios_espectrales[1] = radio_espectral(xk_1, x, x0, 1)
    xk = algoritmo(matriz, b, xk_1)
    errores[2] = rtol(xk, xk_1)
    k = 2
    while True:
        radios_espectrales[k] = radio_espectral(xk, x, x0, k)
        xk = algoritmo(matriz, b, xk_1)
        errores[k] = rtol(xk, xk_1)
        if rtol(xk, xk_1) < 0.001:
            break
        xk_1 = xk
        k += 1
    print xk # resulución
    return {"errores": errores, "radios_espectrales": radios_espectrales}

def test(matriz, algoritmo):
    """ """
    x = [4, 4, 2]*(matriz[0].size/3)
    b = np.dot(matriz, np.array(x, dtype=np.float64))
    # print b # descomentar para ver
    # los vectores b
    x = resolver(algoritmo, matriz, b, x)
    return x

def tocsv(resultados):
    """Exporta los datos a un archivo CSV"""
    with open('tp.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for r in sorted(resultados.keys()):
            writer.writerow(["A" + str(r),])
            for tipo in sorted(resultados[r].keys()):
                writer.writerow(["con " + tipo,])
                writer.writerow(["k", "errores", "radios espectrales"])
                for k in sorted(resultados[r][tipo]["errores"].keys()):
                    writer.writerow([
                        str(k),
                        resultados[r][tipo]["errores"][k],  
                        resultados[r][tipo]["radios_espectrales"][k]])

if __name__ == '__main__':
    resultados = {}
    for i in (6, 18, 24, 30):
        resultados[i] = {}

        print "matriz A" + str(i)
        matriz = crear_matriz(i)
        
        print "Con jacobi"
        resultados[i] = {"jacobi": test(matriz, jacobi_paso)}

        print "Con GS"
        resultados[i]["gs"] = test(matriz, gs_paso)

        print "-"
    tocsv(resultados)
