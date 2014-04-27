
import numpy as np

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
# Se representa como una matriz de 1x3

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
            #print j
            restar -= 1*x.item((j+i)%fila.size)*fila.item((j+i)%fila.size)

        #restar = -x.item((1+i)%3)*fila.item((1+i)%3) - x.item((2+i)%3)*fila.item((2+i)%3)
        temp = b.item(i) + restar
        x_temp[i] = temp/float(fila.item(i))
        #print x_temp

    return x_temp

def jacobi(matriz, b):

    x = jacobi_paso(matriz, b)
    for i in xrange(100):
        x = jacobi_paso(matriz, b, x)
    return x

def test(matriz):
    b = np.dot(matriz, np.array([4, 4, 2]*(matriz[0].size/3) , dtype=np.float64))
    x = jacobi(matriz, b)
    return x

def crear_matriz(n):
    matriz = np.concatenate((B,E), axis=1)
    matriz = np.concatenate((matriz, np.zeros((3,n-6))), axis=1)
    # print matriz
    # agrego lo del medio
    if n > 6:
        cosnt_medio = np.concatenate((E,C,E), axis=1)
        #print medio
        for i in xrange(0, n-6, 3):

            medio = np.concatenate((np.zeros((3,i)), cosnt_medio), axis=1) if i > 0 else cosnt_medio
            # print i, "i"
            medio = np.concatenate((medio, np.zeros((3,n -(9+i)))), axis=1) if n -(9+i) > 0 else medio
            # print n -(9+i), "n -(9+i)"
            # medio = np.concatenate((np.zeros((3,i)), medio, np.zeros((3,n -(9+i)))), axis=1)
            # print matriz, "matriz", matriz.shape
            # print medio, "medio", medio.shape
            matriz = np.concatenate((matriz, medio), axis=0)
    
    fin  = np.concatenate((E,D), axis=1)
    fin = np.concatenate((np.zeros((3,n-6)), fin), axis=1)

    matriz = np.concatenate((matriz, fin), axis=0)

    return matriz

if __name__ == '__main__':
    # print test(B)
    # print test(E)
    # print test(C)
    # print test(D)
    # print test(A6)

    #print crear_matriz(6)
    #print crear_matriz(9)
    #print crear_matriz(12)
    for i in crear_matriz(18):
        print i.tolist()
    print "-"
    #print A6
