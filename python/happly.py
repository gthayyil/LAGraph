import pygraphblas as gb
from operator import itemgetter
def happly(u,x,alpha):
    '''
    The happly function, applies a Householder Reflection
    y= happly(u,x,alpha)
    calculates y = H*x where H = (I - u*u'/alpha)

    H*x =(I - u*u'/alpha)*x=(x-u*(u'u*x)/alpha)

    y,u, and x are vectors of size n and alpha is a scalar.
    The inputs are all modified so y must be a different vector.
    
    Example:
    I=[4.0,5.0,6.0,7.0]
    v=gb.Vector.from_list(I)
    J=[1.0,2.0,3.0,4.0]
    x=gb.Vector.from_list(J)
    A = 2.0
    k=happly(v,x,A)
    print(k)

    This example, simply creates 2 vectors and uses a literal numerical value as alpha to calculate happly.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
    '''
    y =x-u*((u.emult(x)).reduce_float(mon=gb.types.FP32.PLUS_MONOID)/alpha)   
    '''

    s = sum (x) is almost u'*x
    s = ((u'*x) / alpha) ;
    
    assign scalar to x with accum:
       x -= scalar 
        
    does:
        for i = 1:n
            x (i) = x(i) - scalar
    '''
    return y

I=[4.0,5.0,6.0,7.0]
u=gb.Vector.from_list(I)
J=[1.0,2.0,3.0,4.0]
x=gb.Vector.from_list(J)
A = 2.0
gb.options_set(burble=True)
k=happly(u,x,A)
print(k)
