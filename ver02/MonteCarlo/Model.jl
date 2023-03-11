
"""
Electric Interaction
"""
struct K_Model{I<:Integer,F<:Real}
    N::I
    β::F
    mₑ::F
end

function K_Modelly(N::Int64,T::Float64,U::Unit{Float64})
    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    K_Model(N,β,mₑ)
end

function (Problem::K_Model)(φ)
    @unpack N, β, mₑ = Problem
    return H(reshape(φ,3,N),N,mₑ)
end
