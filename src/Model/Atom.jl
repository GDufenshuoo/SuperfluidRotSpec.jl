"""
## Electric Interaction of Atom
`N`: Number of Particle
`B`: Beads of Path
`Z`: Atomic Number
`T`: Temperature/K
`U`: Atomic unit (default)
"""
struct Atom_Model{I<:Integer,F<:Real}
    N::I
    B::I
    β::F
    Z::I
    eꜛ::I
end

function Atom_Model(N::Int64, B::Int64, Z::Int, eꜛ::Int, T::Float64;U::Unit{Float64}=Atomicᵁ)
    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    Atom_Model(N,B,β,Z,eꜛ)
end

function (Problem::Atom_Model)(φ)
    @unpack N, B, β, Z, eꜛ = Problem
    E = 𝑇ᴱ(reshape(φ,3,B,N),N,B,β,eꜛ) + 
        𝑈(reshape(φ,3,B,N),N,B,Z)*β
    return -E
end

"""
# The part to simulate fermions
"""
function 𝑇ᴱ(x,N::Int,B::Int,β::Real,eꜛ::Int)
    T = 0.0
    k = B/β

    for b in 1:B
        T += det(AD(x[:,:,1:eꜛ],eꜛ,B,b,k)) +
            det(AD(x[:,:,eꜛ+1:N],N-eꜛ,B,b,k))
    end
    return -log(T^2)
end

function 𝑈(x,N::Int,B::Int,Z::Int)
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
