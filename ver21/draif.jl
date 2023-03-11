using DynamicHMC, ForwardDiff, LogDensityProblemsAD, Random
using UnPack
using TransformedLogDensities,TransformVariables

using LinearAlgebra
using Zygote

include("Freemove.jl")

struct Unit{F<:Real}
    cᵁ::F
    ħ::F
    eᶜ::F
    mₑ::F
    a₀::F
    sᵁ::F
    Eᵁ::F

    Eᵁₖ::F
end

Atomicᵁ = Unit(1.,1.,1.,1.,1.,1.,1.,3.1668105084779793e-6)

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
    βE = 𝑇ᴱ(reshape(φ,3,B,N),N,B,β)
    return -βE
end

Problem = K_Modelly(2, 12, 10000.0, Atomicᵁ)

ℓ_dims = 3*Problem.N*Problem.B
T = as(Array, ℓ_dims);
ℓ = TransformedLogDensity(T, Problem);  
∇ℓ = ADgradient(:Zygote, ℓ);

rng = Random.GLOBAL_RNG

wu = DynamicHMC.default_warmup_stages()

function extract_initialization(state)
    (; Q, κ, ϵ) = state.final_warmup_state
    (; q = Q.q, κ, ϵ)
end

state1 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[1:1])
state2 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[2:2],
                                     initialization = extract_initialization(state1))
# state3 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[3:3],
#                                      initialization = extract_initialization(state2))
# state4 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[4:4],
#                                      initialization = extract_initialization(state3))
# state5 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[5:5],
#                                      initialization = extract_initialization(state4))
                                
# just keep doing this, and run the last stage with as many samples as you needr