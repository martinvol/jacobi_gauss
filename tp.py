
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

VECTOR_ORIGEN = np.array([0,0,0], dtype=np.float64,)
# Se representa como una matriz de 1x3

def jacobi(matrix, b, x=VECTOR_ORIGEN):
    """Matriz es una matriz de 3x3,
    b es un vector 3 y x un vector semilla de 3 """

    x_temp = np.array(x, dtype=np.float64)

    for i, fila in enumerate(matrix):
        #print i
        temp = b.item(i) - x.item((1+i)%3)*fila.item((1+i)%3) - x.item((2+i)%3)*fila.item((2+i)%3)
        # print x_temp[i]
        a = temp/float(fila.item(i))
        # print a
        x_temp[i] = a
        # print x_temp
        # print x_temp[i], "resultado"

    return x_temp


if __name__ == '__main__':
    x = VECTOR_ORIGEN
    for i in range(100):
        x = jacobi(B, np.array([4,4,2], dtype=np.float64), x)
        print x
    print np.dot(B, x)