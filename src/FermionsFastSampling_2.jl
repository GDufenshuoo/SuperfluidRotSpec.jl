
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
@benchmark factorize([]
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


function FFs_Gauss(A)
    # ith loop j
    for i in 1:N-1, j in i+1:N
        A[:,j] .-= A[:,i].*(A[i,j]/A[i,i])
    end
    for i in N:-1:2, j in i-1:-1:1
        A[:,j] .-= A[:,i].*(A[i,j]/A[i,i])
    end
    return A
end

function 𝑝(A,B)
    sum(abs2,A.-B)
end
