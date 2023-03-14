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
    E = 𝑇ᴱ_Atom(reshape(φ,3,B,N),N,B,β,eꜛ) + 
        𝑈_Atom(reshape(φ,3,B,N),N,B,Z)
    return -E
end

"""
# The part to simulate fermions
"""
function 𝑇ᴱ_Atom(x,N::Int,B::Int,β::Real,eꜛ::Int)
    T = 0.0
    k = -0.5*B/β

    for b in 1:B
        T += AD(x[:,:,1:eꜛ],eꜛ,B,b,k) + 
            AD(x[:,:,eꜛ+1:N],N-eꜛ,B,b,k)
    end
    return -log(abs(T))*2/β
end

function 𝑈_Atom(x,N::Int,B::Int,Z::Int)
    U = 0.0
    @floop for i in 1:N
        for b in 1:B
            @reduce U += Z/(norm(x[:,b,i])+1e-10)
    end end
    @floop for i in 2:N
        for j in 1:i, b in 1:B
            @reduce U -= 1/(norm(x[:,b,i].-x[:,b,j])+1e-10)
    end end
    return U/B
end
