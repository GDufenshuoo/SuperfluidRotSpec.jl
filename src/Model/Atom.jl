# Need fix

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

function Atom_Model(N::Int, B::Int, Z::Int, eꜛ::Int, T::Float64; U::Unit{Float64}=Atomicᵁ)
    @unpack mₑ, ħ, Eᵁₖ = U
    Atom_Model(N,B,T,Z,eꜛ)
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
    return -log(abs(T))*2β
end

function 𝑈_Atom(x,N::Int,B::Int,Z::Int)
    U1 = 0.0
    U2 = 0.0
    for i in 1:N
        for b in 1:B
            U1 -= 1/norm(x[:,b,i])
    end end
    for i in 2:N
        for j in 1:i-1, b in 1:B
            U2 += 1/norm(x[:,b,i].-x[:,b,j])
    end end
    return (Z*U1+U2)/B
end
