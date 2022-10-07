from happly import happly
from norm2 import norm2
import pygraphblas as gb
def hmhx(M,u,x,alpha):
    '''
    hmhx = a function used to compute y = H*M*H*x = (M-u*x'-x*u)*x where x and y have the same size n 

    Example:
    I=[4.0,5.0,6.0,7.0]
    v=gb.Vector.from_list(I)
    x=gb.Vector.from_list(I)
    A = 2.0
    I=[3,3,3]
    J=[1,2,3]
    K=[4,5,6]
    M=gb.Matrix.from_lists(I,J,K)
    M=M+M.transpose()
    mdiag=M.diag()
    indiag=1/mdiag.cast(gb.FP64)
    print(hmhx(indiag,v,x,A))
    
    This example starts off by creating 2 vectors and a matrix. 
    Then makes sure the matrix is symmetric, and  gets the inverse diagonal of the matrix. 
    Lastly, the example prints, hmhx called using the modified matrix, and vectors.


    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
    
    '''
    print("---------------------------------------------------------------------------------------------------------------")
    print("ENTERING HMHX")
    print("---------------------------------------------------------------------------------------------------------------")
    #Calculates the happly with u,x, and alpha
    y=happly(u,x,alpha)
    #Matrix vector multiply using M and y
    g=M.mxv(y)
    #Calculates happly with u,g, and alpha
    z=happly(u,g,alpha)
    #Though z[0] is bound to be a very low value, manually sets it to 0, to prevent issues. 
    z[0]=0
    #Returns z
    print("---------------------------------------------------------------------------------------------------------------")
    print("ABOUT TO EXIT HMHX")
    print("---------------------------------------------------------------------------------------------------------------")
    return z


