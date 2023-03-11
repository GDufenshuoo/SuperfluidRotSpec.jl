
"""
Electric Interaction
"""
struct K_Model{I<:Integer,F<:Real}
    N::I
    B::I
    β::F
end

function K_Modelly(N::Int64, B::Int64, T::Float64, U::Unit{Float64})
    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    K_Model(N,B,β)
end

function (Problem::K_Model)(φ)
    @unpack N, B, β = Problem
    return 𝑈(reshape(φ,3,B,N)) + 𝑇ᴱ(reshape(φ,3,B,N))
end
