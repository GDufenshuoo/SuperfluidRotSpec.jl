"""
## Electric Interaction of Atom
`N`: Number of Particle
`B`: Beads of Path
`Z`: Atomic Number
"""
struct Atom_Model{I<:Integer,F<:Real}
    N::I
    B::I
    β::F
    Z::I
    eꜛ::I
end


"""
# Set the Model
`N`: Number of Particle
`B`: Beads of Path
`T`: Temperature/K

`U`: Atomic unit (default)
"""
function Atom_Model(N::Int64, B::Int64, T::Float64, Z::Int, eꜛ::Int;U::Unit{Float64}=Atomicᵁ)
    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    K_Model(N,B,β,Z,eꜛ)
end

function (Problem::Atom_Model)(φ)
    @unpack N, B, β = Problem
    E = 𝑇ᴱ(reshape(φ,3,B,N),Problem::Atom_Model) + 
        𝑈(reshape(φ,3,B,N),Problem::Atom_Model)
    return -E
end

"""
# The part to simulate fermions
"""
function 𝑇ᴱ(x,Problem::Atom_Model)
    @unpack N,B,β,eꜛ = Problem
    T = 0.0
    k = B/β

    for b in 1:B
        T += det(AD(x[:,:,1:eꜛ],N,B,b,k)) +
            det(AD(x[:,:,eꜛ+1:N],N,B,b,k))
    end
    return -log(T^2)/β
end

function AD(x,N,B,b,k)
    A = Zygote.Buffer(zeros(),N, N)
    L = (b == 1 ? B : b-1)
    for i in 1:N, j in 1:N
        A[i,j] = exp(-0.5/k*𝑝(x[:,L,i],x[:,b,j]))
    end
    return copy(A)
end

function 𝑈(x,Problem::Atom_Model)
    @unpack N,B,Z = Problem
    U = 0.0
    for i in 1:N
        for b in 1:B
            U += Z/(norm(x[:,b,i])+1e-10)
    end end
    for i in 2:N
        for j in 1:i, b in 1:B
            U -= 1/(norm(x[:,b,i].-x[:,b,j])+1e-10)
    end end
    return U/B
end
