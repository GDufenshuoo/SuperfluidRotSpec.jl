function extract_initialization(state)
    (; Q, κ, ϵ) = state.final_warmup_state
    (; q = Q.q, κ, ϵ)
end

"""
## run HMC to solve the problem
"""
function runHMC(Problem,warmup::Int,Num::Int)
    @unpack N,B = Problem
    println(show(Problem),"\n Begin to HMC")
    ℓ_dims = 3*N*B
    T = as(Array, ℓ_dims);
    ℓ = TransformedLogDensity(T, Problem);  
    ∇ℓ = ADgradient(:Zygote, ℓ);
    rng = Random.GLOBAL_RNG
    wu = DynamicHMC.default_warmup_stages()

    state, warmup_stage = HMC_warmup(rng,∇ℓ,warmup,wu)
    
    return DynamicHMC.mcmc_with_warmup(rng, ∇ℓ, Num;
            warmup_stages = wu[warmup_stage:warmup_stage],
            initialization = extract_initialization(state))
end

"""
## run HMC to solve the problem
"""
function runHMC(Problem,Num::Int)
    @unpack N,B = Problem
    println(show(Problem),"\n Begin to HMC")
    ℓ_dims = 3*N*B
    T = as(Array, ℓ_dims);
    ℓ = TransformedLogDensity(T, Problem);  
    ∇ℓ = ADgradient(:Zygote, ℓ);
    rng = Random.GLOBAL_RNG
    return DynamicHMC.mcmc_KEEP_warmup(rng, ∇ℓ, Num;)
end

