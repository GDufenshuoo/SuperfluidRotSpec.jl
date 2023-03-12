module SuperfluidRotSpec

export 𝑇ᴱ
export K_Modelly, Atomicᵁ, Unit
export PI

using DynamicHMC, ForwardDiff, LogDensityProblemsAD, Random
using UnPack
using TransformedLogDensities,TransformVariables
using FLoops
using LinearAlgebra
using Zygote

include("2019Bosonic.jl")
include("PE.jl")

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

const Atomicᵁ = Unit(1.,1.,1.,1.,1.,1.,1.,3.1668105084779793e-6)

"""
Electric Interaction
"""
struct K_Model{I<:Integer,F<:Real}
    N::I
    B::I
    β::F
end

"""
# Set the Model
`N`: Number of Particle
`B`: Beads of Path
`T`: Temperature/K

`U`: Atomic unit default
"""
function K_Modelly(N::Int64, B::Int64, T::Float64;U::Unit{Float64}=Atomicᵁ)
    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    K_Model(N,B,β)
end

function (Problem::K_Model)(φ)
    @unpack N, B, β = Problem
    E = 𝑇ᴱ(reshape(φ,3,B,N),N,B,β)/β + 𝑈(reshape(φ,3,B,N),N,B)
    return -E
end

function extract_initialization(state)
    (; Q, κ, ϵ) = state.final_warmup_state
    (; q = Q.q, κ, ϵ)
end

"""

"""
function PI(Problem::K_Model,warmup::Int,Num::Int)
    @unpack N,B = Problem
    ℓ_dims = 3*N*B
    T = as(Array, ℓ_dims);
    ℓ = TransformedLogDensity(T, Problem);  
    ∇ℓ = ADgradient(:Zygote, ℓ);
    rng = Random.GLOBAL_RNG
    state,warmup_stage = HMC_warmup(rng,∇ℓ,warmup)
    return DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, Num;
            warmup_stages = wu[warmup_stage:warmup_stage],
            initialization = extract_initialization(state))
end

function HMC_warmup(rng,∇ℓ,warmup::Int)
    # just keep doing this, and run the last stage with as many samples as you need
    if warmup > 5
        println("warmup_stages should less then 6: Auto-set warmup = 5")
        warmup = 5
    end
    wu = DynamicHMC.default_warmup_stages()
    state = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[1:1])
    for i in 1:warmup-1
        state = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[i+1:i+1],
                                            initialization = extract_initialization(state))
    end
    return stage, warmup+1
end

end
