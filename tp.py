
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

def jacobi_paso(matrix, b, x=VECTOR_ORIGEN):
    """Matriz es una matriz de 3x3,
    b es un vector 3 y x un vector semilla de 3 """


    x_temp = np.array(x, dtype=np.float64)

    for i, fila in enumerate(matrix):

        # print "x%d = b%d - a%d%dX%d - a%d%dX%d /%d%d" % (i, i, i, (1+i)%3, (1+i)%3, i, (2+i)%3, (2+i)%3, i, i)
        # print "x%d = %d - %d*%d - %d*%d /%d" % (i, b.item(i), x.item((1+i)%3), fila.item((1+i)%3), x.item((2+i)%3), fila.item((2+i)%3), float(fila.item(i)))
        temp = b.item(i) - x.item((1+i)%3)*fila.item((1+i)%3) - x.item((2+i)%3)*fila.item((2+i)%3)
        x_temp[i] = temp/float(fila.item(i))
        #print x_temp

    return x_temp

def jacobi(matriz, b):

    x = VECTOR_ORIGEN
    for i in range(100):
        x = jacobi_paso(matriz, b, x)
    return x

def test(matriz):
    b = np.dot(matriz, np.array([4, 4, 2], dtype=np.float64))
    x = jacobi(matriz, b)
    return x#np.dot(matriz, x)

if __name__ == '__main__':
    print test(B)
    print test(E)
    print test(C)
    print test(D)
