"""
## run HMC to solve the problem
"""
function runHMC(ℓ::ClassicRotor;
    lP_samples = 10_000, lP_adapts = 5_000, initθ)
    @unpack N = ℓ
    println("Classic Rotor \n Begin to HMC")
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
    lP, lP_stats = sample(rng, hamiltonian, proposal, initθ, lP_samples, adaptor, lP_adapts;progress=true)
    return lP[1:lP_adapts], lP[end-lP_adapts+1:end], lP_stats
end

"""
## run HMC to solve the problem
"""
function runHMC(ℓ::SuperfluidRotor;
    lP_samples = 10_000, lP_adapts = 5_000, initθ)
    @unpack N, B = ℓ
    println("Superfluid Rotor \n Begin to HMC")
    rng = Random.GLOBAL_RNG
    T = as(Array, 3*N*B)
    ℓ = TransformedLogDensity(T, ℓ)
    metric = DiagEuclideanMetric(3N*B)
    hamiltonian = Hamiltonian(metric, ℓ, ReverseDiff)
    initial_ϵ = find_good_stepsize(hamiltonian, initθ)
    integrator = Leapfrog(initial_ϵ)
    proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
    adaptor = StanHMCAdaptor(rng, MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
    lP, lP_stats = sample(hamiltonian, proposal, initθ, lP_samples, adaptor, lP_adapts;progress=true)
    return lP[1:lP_adapts], lP[end-lP_adapts:end], lP_stats
end