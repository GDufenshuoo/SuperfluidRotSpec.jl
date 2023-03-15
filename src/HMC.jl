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

"""
It Can't be more usefull (although not a very good code)
"""
function HMC_warmup(rng,∇ℓ,warmup::Int,wu)
    # just keep doing this, and run the last stage with as many samples as you need
    if warmup > 5
        println("warmup_stages should less then 6: Auto-set warmup = 5")
        warmup = 5
    end
    i = 1
    println(">>> Stage $i ...")

    state = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[1:1])
    println("Finished \n")
    if warmup == i return state,i+1 else i+=1 end
    
    println(">>> Stage $i ...")
    state = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[i:i],
                                        initialization = extract_initialization(state))
    println(">>> Stage $i Finished \n")
    if warmup == i return state,i+1 else i+=1 end

    println(">>> Stage $i ...")
    state = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[i:i],
                                        initialization = extract_initialization(state))
    println(">>> Stage $i Finished \n")
    if warmup == i return state,i+1 else i+=1 end

    println(">>> Stage $i ...")
    state = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[i:i],
                                        initialization = extract_initialization(state))
    println(">>> Stage $i Finished \n")
    if warmup == i return state,i+1 else i+=1 end

    println(">>> Stage $i ...")
    state = DynamicHMC.mcmc_keep_warmup(rng, ∇ℓ, 0; warmup_stages = wu[i:i],
                                        initialization = extract_initialization(state))
    println(">>> Stage $i Finished \n")

    return state,i+1
end

