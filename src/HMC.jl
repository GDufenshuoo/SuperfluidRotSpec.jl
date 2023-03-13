function extract_initialization(state)
    (; Q, κ, ϵ) = state.final_warmup_state
    (; q = Q.q, κ, ϵ)
end

"""
## run HMC to solve the problem
"""
function runHMC(Problem,warmup::Int,Num::Int)
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

