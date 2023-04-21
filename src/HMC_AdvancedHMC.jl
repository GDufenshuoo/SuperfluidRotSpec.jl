"""
## run HMC to solve the problem
"""
function runHMC(ℓ::ClassicRotor;
    lP_samples = 10_000, lP_adapts = 5_000, initθ, showpro=false)
    @unpack N = ℓ
    println(">"^10*"Classic Rotor \n Begin to HMC")
    rng = Random.GLOBAL_RNG

    T = as(Array, 3*N)
    ℓ = TransformedLogDensity(T, ℓ)
    metric = DiagEuclideanMetric(3N)
    hamiltonian = Hamiltonian(metric, ℓ, ReverseDiff)
    initial_ϵ = find_good_stepsize(hamiltonian, initθ)
    integrator = Leapfrog(initial_ϵ)
    proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
    # adaptor = StanHMCAdaptor(StepSizeAdaptor(0.8, integrator))
    lP, lP_stats = sample(rng, hamiltonian, proposal, initθ, lP_samples, adaptor, lP_adapts;progress=showpro)
    return lP[1:lP_adapts], lP[end-lP_adapts+1:end], lP_stats
end

"""
## run HMC to solve the problem
"""
function runHMC(ℓ::SuperfluidFixRotor;
    lP_samples = 10_000, lP_adapts = 5_000, initθ, showpro=false)
    @unpack N, B = ℓ
    println(">"^10*"Superfluid fixed-Rotor \n Begin to HMC")
    rng = Random.GLOBAL_RNG
    T = as(Array, 3*N*B)
    ℓ = TransformedLogDensity(T, ℓ)
    metric = DiagEuclideanMetric(3*N*B)
    hamiltonian = Hamiltonian(metric, ℓ, ReverseDiff)
    initial_ϵ = find_good_stepsize(hamiltonian, initθ)
    integrator = Leapfrog(initial_ϵ)
    proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
    lP, lP_stats = sample(hamiltonian, proposal, initθ, lP_samples, adaptor, lP_adapts;progress=showpro)
    return lP[1:lP_adapts], lP[end-lP_adapts+1:end], lP_stats
end

"""
## run HMC to solve the problem
"""
function runHMC(ℓ::SuperfluidRotor;
    lP_samples = 10_000, lP_adapts = 5_000, initθ, showpro=false)
    @unpack N, B, rRB = ℓ
    RB = fld(B,rRB)
    println(">"^10*"Superfluid Rotor \n Begin to HMC")
    rng = Random.GLOBAL_RNG
    T = as(Array, 3*N*B+5*RB)
    ℓ = TransformedLogDensity(T, ℓ)
    metric = DiagEuclideanMetric(3*N*B+5*RB)
    hamiltonian = Hamiltonian(metric, ℓ, ReverseDiff)
    initial_ϵ = find_good_stepsize(hamiltonian, initθ)
    integrator = Leapfrog(initial_ϵ)
    proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
    lP, lP_stats = sample(hamiltonian, proposal, initθ, lP_samples, adaptor, lP_adapts;progress=showpro)
    return lP[1:lP_adapts], lP[end-lP_adapts+1:end], lP_stats
end