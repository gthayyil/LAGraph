import pygraphblas as gb
from operator import itemgetter
import time
#gb.options_set(burble=True)
def laplacian(G):
    '''
    Laplacian: This function computes the laplacian of a Matrix (This matrix has to be symmetric (a normal matrix, when added with its transpose becomes symmetric)), and returns a tuple with 
    both the laplacian of the matrix, and the infinity norm of the matrix.

    Example:
    J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
    h=J[1]+J[1].transpose()
    L=laplacian(h)
    print(L)

    The first line is simply the way to get a matrix from the SuiteSparse collection in python using the
    id of a matrix. It requires the import of "itemgetter" from operator, and the prior hand installation
    of ssgetpy.

    The second line simply makes sure that the matrix is symmetric by adding it with its transpose.

    The third line calls this laplacian function on the matrix and returns a tuple with the laplacian and the inf norm.

    Authors: Georgy Thayyil, Tim Davis
    Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm
    '''
    print("-----------------------------------------------------------------------------------------------")
    print("INSIDE LAPLACIAN FUNCTION")
    print("-----------------------------------------------------------------------------------------------")

    #The offdiag function is used to create a matrix that doesnt have the diagonal of the matrix G.
    A=G.offdiag()  

    #Created a dense vector filled with 0s, with the size of the number of columns of matrix A.
    d=gb.Vector.dense(gb.types.FP32,A.ncols) 

    #Uses a matrix vector multiply with the PLUS_PAIR semiring 
    # d += A*d  
    A.mxv(d,accum=gb.types.FP32.PLUS,semiring=gb.types.FP32.PLUS_PAIR,out=d)

    #Creates a sparse Matrix with same dimensions as A, and asssigns -1 with the A as a Mask.
    S=gb.Matrix.sparse(gb.types.FP32,A.nrows,A.ncols)
    S.assign_scalar(-1,mask=A,desc=gb.descriptor.S)
    
    #Finds all the places where 0 appear in the passed in matrix, and uses that as a mask to replace those places with 1.   
    k=d.select('==',0)
    d.assign_scalar(1,mask=k,desc=gb.descriptor.S)

    #Find the inf norm using the vector d and MAX_MONOID
    inform=d.reduce_float(d.type.MAX_MONOID)*2

    #Using a built in function to create a diagonal matrix from a vector.
    Dmatrix=gb.Matrix.from_diag(d)     

    #Calculating the Laplacian by adding the Dmatrix with S. 
    L=Dmatrix+S

    print("-----------------------------------------------------------------------------------------------")
    print("ABOUT TO EXIT LAPLACIAN FUNCTION")
    print("-----------------------------------------------------------------------------------------------")
    #Returns the laplacian and the inf norm as a tuple
    return [L,inform]
J=sorted(list(gb.Matrix.ssget(2760)), key=itemgetter(0))[0]
h=J[1]+J[1].transpose()
L=laplacian(h)
print(L)

