using AdvancedHMC
using ForwardDiff
using Zygote
using LogDensityProblems
using LinearAlgebra

include("Bosonic.jl")


P = ones(3,B,N)

A = rand(3,N)

for i in 1:N,j in 1:B
    P[:,j,i] = A[:,i] .+ 0.01*randn(3)
end

E_guess = 2000

# Define the target distribution using the `LogDensityProblem` interface
struct LogTargetDensity
    dim::Int
end
LogDensityProblems.logdensity(p::LogTargetDensity, θ) = log.(WᴮO(reshape(θ,3,B,N)))   # standard multivariate normal
LogDensityProblems.dimension(p::LogTargetDensity) = p.dim
LogDensityProblems.capabilities(::Type{LogTargetDensity}) = LogDensityProblems.LogDensityOrder{0}()

# Choose parameter dimensionality and initial parameter value
D = 3*N*B; initial_θ = P[:]
ℓπ = LogTargetDensity(D)

<<<<<<< Updated upstream
# Set the number of samples to draw and warmup iterations
n_samples, n_adapts = 100, 50

# Define a Hamiltonian system
metric = DiagEuclideanMetric(D)
hamiltonian = Hamiltonian(metric, ℓπ, Zygote)

# Define a leapfrog solver, with initial step size chosen heuristically
initial_ϵ = find_good_stepsize(hamiltonian, initial_θ)
integrator = Leapfrog(initial_ϵ)
TemperedLeapfrog(initial_ϵ, 10.0)

# Define an HMC sampler, with the following components
#   - multinomial sampling scheme,
#   - generalised No-U-Turn criteria, and
#   - windowed adaption for step-size and diagonal mass matrix
# proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
proposal = NUTS{SliceTS,GeneralisedNoUTurn}(integrator)
# proposal = HMCDA(integrator, 1.0)
adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))

# Run the sampler to draw samples from the specified Gaussian, where
#   - `samples` will store the samples
#   - `stats` will store diagnostic statistics for each sample
samples, stats = AdvancedHMC.sample(hamiltonian, proposal, initial_θ, n_samples, adaptor, n_adapts; progress=true)


# using Plots

# plot(samples[1000][:])
=======
>>>>>>> Stashed changes
