
import numpy as np
from numpy.core import Inf

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
    """Matriz es una matriz de 3x3,
    b es un vector 3 y x un vector semilla de 3 """
    if x is None:
        x = np.zeros(matrix[0].size)

    x_temp = np.array(x, dtype=np.float64)

    for i, fila in enumerate(matrix):
        # print "x%d = b%d - a%d%dX%d - a%d%dX%d /%d%d" % (i, i, i, (1+i)%3, (1+i)%3, i, (2+i)%3, (2+i)%3, i, i)
        # print "x%d = %d - %d*%d - %d*%d /%d" % (i, b.item(i), x.item((1+i)%3), fila.item((1+i)%3), x.item((2+i)%3), fila.item((2+i)%3), float(fila.item(i)))
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
        # print x_temp
        # raw_input()

    return x_temp

def gs_paso(matrix, b, x=None):
    """Matriz es una matriz de 3x3,
    b es un vector 3 y x un vector semilla de 3 """
    if x is None:
        x = np.zeros(matrix[0].size)

    x_temp = np.array(x, dtype=np.float64)

    for i, fila in enumerate(matrix):
        # print "x%d = b%d - a%d%dX%d - a%d%dX%d /%d%d" % (i, i, i, (1+i)%3, (1+i)%3, i, (2+i)%3, (2+i)%3, i, i)
        # print "x%d = %d - %d*%d - %d*%d /%d" % (i, b.item(i), x.item((1+i)%3), fila.item((1+i)%3), x.item((2+i)%3), fila.item((2+i)%3), float(fila.item(i)))
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
        # print x_temp
        # raw_input()

    return x_temp


def rtol(xk, xk_1):

    print np.linalg.norm(xk - xk_1 , ord=Inf)/np.linalg.norm(xk, ord=Inf)
    return np.linalg.norm(xk - xk_1 , ord=Inf)/np.linalg.norm(xk, ord=Inf)

def resolver(algoritmo, matriz, b):
    xk_1 = jacobi_paso(matriz, b)
    xk = jacobi_paso(matriz, b, xk_1)
    while True:
        xk = algoritmo(matriz, b, xk_1)
        if rtol(xk, xk_1) < 0.001:
            break
        xk_1 = xk
    return xk

def test(matriz, algoritmo):
    b = np.dot(matriz, np.array([4, 4, 2]*(matriz[0].size/3) , dtype=np.float64))
    x = resolver(algoritmo, matriz, b)
    return x

if __name__ == '__main__':
    matriz = crear_matriz(30)
    # print matriz
    print test(matriz, jacobi_paso)
    print test(matriz, gs_paso)
