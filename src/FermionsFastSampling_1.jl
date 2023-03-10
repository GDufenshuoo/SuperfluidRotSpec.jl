
"""

"""
function Samplingᵒ(Pos_pre,Pos_now;Mcond=1e5)
    Pos_new = copy(Pos_now)
    E = Eᵒ(Pos_pre,Pos_now)
    move = randn(3).*Mᵥ

    for i in 1:Mcond
        perm = 
        for n in 1:N
            move = randn(3).*Mᵥ
            Pos_new[:,n] .+= move
            if Accᵒ(Pos_pre, Pos_now, Pos_new,perm)#; Fermions = true)
                print("O",move)
            else
                Pos_new[:,n] .-= move
    end end end
    return Pos_new
end

𝑘 = 1
N = 2

Pos_pre = 10*rand(3,2)
Pos_now = 10*rand(3,2)

Mᵥ = 1
a = Samplingᵒ(Pos_pre,Pos_now;Mcond=1e3).-Pos_now

A = zeros(N,N)

FFs(randn(N,N))

B1 = 1.0
C1 = FFs(randn(N,N))
C2 = FFs(randn(N,N))
for i in 1:N
    B1 *= C1[i]/C2[i]
end
B1^2

N = 200

B = FFs(rand(3,N),rand(3,N))
# const 
𝑘 = 1

using BenchmarkTools
using LinearAlgebra
@benchmark factorize(B)
@benchmark FFs_Gauss(B)

factorize(B)
FFs_Gauss(B)


function FFs(px,x)
    A = zeros(N,N)
    for i in 1:N, j in 1:N
        A[i,j] = exp(-𝑘*𝑝(px[:,i],x[:,j]))
    end
    return A
end

"""
import numpy as np

def calculate_normal_vector(U, k):
    ""
    Calculate the normal vector for the k-th row of the matrix U, given the
    first k-1 rows have been transformed to reduced row echelon form using
    the iterative Gaussian elimination method.
    
    Args:
    U (numpy array): An Nxk matrix whose first k-1 rows have been transformed
        to reduced row echelon form.
    k (int): The index of the row to be transformed next.
    
    Returns:
    numpy array: The normal vector for the k-th row of U.
    ""
    # Copy the first k-1 rows of U to a new matrix and add a row of zeros for the
    # k-th row (which has not been transformed yet).
    U_temp = np.zeros((k, k))
    U_temp[:-1,:] = U[:k-1,:]
    
    # Perform the next step of Gaussian elimination to transform the k-th row.
    for i in range(k-1):
        pivot = U_temp[i,i]
        U_temp[i+1:,i:] -= U_temp[i,i:]*U_temp[i+1:,i]/pivot
        U_temp[i+1:,i] = 0
    U_temp[-1,-1] = 1
    
    # Extract the normal vector from the last row of the transformed matrix.
    normal_vector = U_temp[-1,:k-1]
    
    return normal_vector


def calculate_determinant(U):
    ""
    Calculate the determinant of the matrix U using the volume interpretation.
    
    Args:
    U (numpy array): An NxN matrix.
    
    Returns:
    float: The determinant of U.
    ""
    # Calculate the normal vectors for each row of U using the iterative
    # Gaussian elimination method.
    normals = []
    for k in range(1, U.shape[0]+1):
        normals.append(calculate_normal_vector(U, k))
    
    # Calculate the height of the parallelepiped spanned by the row vectors
    # using the dot product of the last row vector with each normal vector.
    height = np.dot(U[-1,:], normals[-1])
    
    # Calculate the determinant as the absolute value of the volume.
    volume = np.abs(np.linalg.det(normals))
    determinant = volume * height
    
    return determinant
"""


function FFs_Gauss(A)
    # ith loop j
    for i in 1:N-1, j in i+1:N
        A[:,j] .-= A[:,i].*(A[i,j]/A[i,i])
    end
    for i in N:-1:2, j in i-1:-1:1
        A[:,j] .-= A[:,i].*(A[i,j]/A[i,i])
    end
    return A#[A[i,i] for i in 1:N]
end

function 𝑝(A,B)
    sum(abs2,A.-B)
end
