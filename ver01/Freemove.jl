
"""
### Monte Carlo 解自由粒子
"""
function T(x)
    T = 0.0
    for b in 2:B
        for i in 1:N
            T += 𝑝(x[b-1,i],x[b,i])
    end end
    T += 𝑝(x[begin,i],x[end,i])
    return T
end

"""
### get the 𝑝
"""
function 𝑝(a::Vector,b::Vector)
    sum(abs2,a.-b)
end

function Tᵒ(x,b::Int)
    T = 0.0
    L = (b == 1 ? B : b-1)
    for i in 1:N
        T += 𝑝(x[L,i],x[b,i])
    end
    return T
end

"""
### Moving Forward
"""
function Tᴬ(x)
    A = zeros(B,N)
    L = (b == 1 ? B : b-1)
    for b in 1:B
        for i in 1:N
            A[b,i] = 𝑝(x[L,i],x[b,i])
    end end
    return A
end

"""
### A list of moving particles needed
"""
function δTᵒ!(x,RM,b::Int,io::Vector)
    T = 0.0
    L = (b == 1 ? B : b-1)
    for i in eachindex(io)
        n = io[i]
        T += 𝑝(x[L,n],x[b,n].+RM[i])
    end
    return T
end

using LinearAlgebra

"""
### 用 Slater 行列式来解自由费米粒子
"""
function T(x)
    T = 0.0
    A = zeros(N,N)
    for b in 2:B
        for i in 1:N, j in 1:N
            A[i,j] = 𝑝(x[b-1,i],x[b,j])
        end
        T += det(A)
    end
    for i in 1:N, j in 1:N
        A[i,j] = 𝑝(x[begin,i],x[end,j])
    end
    T += det(A)
    return T^2
end

"""
### Moving Forward
"""
function Tᵒ(x,b::Int)
    T = 0.0
    A = zeros(N,N)
    L = (b == 1 ? B : b-1)
    for i in 1:N, j in 1:N
        A[i,j] = 𝑝(x[L,i],x[b,j])
    end
    T += det(A)
    return T^2
end

"""
### Moving Forward
"""
function Tᴬ(x)
    A = zeros(N,N,B)
    for b in 1:B
        L = (b == 1 ? B : b-1)
        for i in 1:N, j in 1:N
            A[i,j,b] = 𝑝(x[L,i],x[b,j])
    end end
    return A
end

"""
### Moving Forward
"""
function Tᴬᵒ(x,b::Int)
    A = zeros(N,N)
    L = (b == 1 ? B : b-1)
    for i in 1:N, j in 1:N
        A[i,j] = 𝑝(x[L,i],x[b,j])
    end
    return A
end

"""
### A list of moving particles needed
"""
function δTᵒ(x,RM,b::Int,io::Vector,C::Matrix)
    A = zeros(N,N)
    L = (b == 1 ? B : b-1)
    x[b,:] += RM
    for i in 1:N, j in 1:N
        if i in io || j in io
            A[i,j] = 𝑝(x[L,i],x[b,j])
        end
    end
    return (det(A+C)-det(C))^2, A+C
end

function T(x::Boses)
    
end
