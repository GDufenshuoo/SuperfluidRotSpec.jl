using DynamicHMC, ForwardDiff, LogDensityProblemsAD, Random
using UnPack
using TransformedLogDensities,TransformVariables

Problem = K_Modelly(2, 240, 1.0, Atomicᵁ)

ℓ_dims = 3*Problem.N
T = as(Array, ℓ_dims);
ℓ = TransformedLogDensity(T, Problem);  
∇ℓ = ADgradient(:ForwardDiff, ℓ);

rng = Random.GLOBAL_RNG

wu = DynamicHMC.default_warmup_stages()

function extract_initialization(state)
    (; Q, κ, ϵ) = state.final_warmup_state
    (; q = Q.q, κ, ϵ)
end

state1 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[1:1])
state2 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[2:2],
                                     initialization = extract_initialization(state1))
state3 = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[3:3],
                                     initialization = extract_initialization(state2))
# just keep doing this, and run the last stage with as many samples as you need