

import pygraphblas as gb
def happly(u,x,alpha):
    y =x-u*((u.emult(x)).reduce_float()/alpha)
    #print(y)
    return y
~               