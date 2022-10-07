import pygraphblas as gb
from math import sqrt
def norm2(v):
    '''
    norm2: This function aims to get the 2 norm of a vetor, and requires an import from python's math module.

    Example:
    I=[4.0,5.0,6.0,7.0]
    v=gb.Vector.from_list(I)
    k=norm2(v)
    print(k)

    Creates a vector, and passes it to the function and gets the norm2 of the function as a result.
    
    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
    '''
    #The passed in vector gets squared
    k=v.apply_second(gb.FP32.POW,2)
    #Now this vector gets reduced to one float (the default is addition)
    z=k.reduce_float(mon=gb.types.FP32.PLUS_MONOID)
    #Now we take the square root of z
    g=sqrt(z)
    #Return the value of norm2 for the passed in vector
    return g

