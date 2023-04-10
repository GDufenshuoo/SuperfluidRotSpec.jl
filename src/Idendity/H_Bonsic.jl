
function Mβ𝐻(X,X̄,τ,λ,m)
    x = reshape(X,2,N,m)
    x̄ = reshape(X̄,2,N,m)
    E = 0.0
    for i in 1:m
        E += β𝐻(x[:,:,i],x̄[:,:,i],τ,λ)
    end
    return E/m
end

"""
H
"""
function β𝐻(x,x̄,τ,λ)
    return β𝑉(x,x̄,τ,λ) + β𝑇ᴱ(x,x̄,τ)
end

function Class_Mβ𝑉(X,λ,m)
    x = reshape(X,2,N,m)
    E = 0.0
    for i in 1:m
        E += 𝑉(x[:,:,i],λ)
    end
    return E/m
end

function Mβ𝑉(X,X̄,τ,λ,m)
    x = reshape(X,2,N,m)
    x̄ = reshape(X̄,2,N,m)
    E = 0.0
    for i in 1:m
        E += β𝑉(x[:,:,i],x̄[:,:,i],τ,λ)
    end
    return E/m
end

function β𝑉(x,x̄,τ,λ)
    return τ*(𝑉(x,λ)+𝑉(x̄,λ))
end

function 𝑉(x,λ)
    return 𝑈(x)/2 + ⋓(x,N,λ)
end

function 𝑈(x)
    U = 0.0
    U += sum(abs2,x)
    return U
end

function ⋓(x,N,λ)
    U = 0.0
    for i in 2:N, j in 1:i-1
        U += λ/sqrt(sum(abs2,x[:,i].-x[:,j]))
    end
    return U
end

"""
# The part to simulate fermions
"""
function β𝑇ᴱ(x,x̄,τ)
    k = -0.5/τ
    T = AD(x,x̄,N,k)
    return -log(T^2)/2
end

"""
Under developing
"""
function AD(x,x̄,N,k)
    E = 0.0
    A = zeros(Real,N, N)
    for i in 1:N
        for j in 1:N
            A[i,j] = exp(k*sum(abs2,x[:,i].-x̄[:,j]))
        end
    end
    E += permanent(A)
    return E
end

function matrix_minor(A::Matrix, row::Int, col::Int)
    minor = [A[i, j] for i in 1:size(A, 1) if i ≠ row for j in 1:size(A, 2) if j ≠ col]
    return reshape(minor, size(A, 1) - 1, size(A, 2) - 1)
end

function permanent(A::Matrix)
    l = size(A, 1)
    if l == 1
        return A[1, 1]
    elseif l == 2
        return A[1, 1] * A[2, 2] + A[1, 2] * A[2, 1]
    else
        det = 0.0
        for col in 1:l
            det += permanent(matrix_minor(A, 1, col))
        end
        return det
    end
end
