# using DynamicHMC
using LinearAlgebra
using Statistics
using AdvancedHMC
using LogDensityProblems
using LogDensityProblemsAD
# using ForwardDiff
using ReverseDiff
using UnPack
include("Idendity/PI.jl")
include("Idendity/2019Bosonic.jl")
# include("Idendity/Fermions.jl")
using JLD2
using Interpolations
include("Units.jl")
include("PE.jl")
include("Model/Model.jl")

# paraH2_paraH2 =  load("OCS.PES","paraH2_paraH2")
# OCS_paraH2 = load("OCS.PES","OCS_paraH2")

# paraH2_paraH2 = set_potention(paraH2_paraH2)
# OCS_paraH2 = set_potention(OCS_paraH2)

N = 2
B = 512
T = 1.0

Problem = SuperfulidRotor(N, B, T, "OCS.PES", "OCS_paraH2", "paraH2_paraH2")
Classic_Problem = Classic_SuperfulidRotor(N, 1, Problem.β, Problem.rotor, Problem.superfulid, Problem.am2An)

struct Classic_SuperfulidRotor{I<:Integer,F<:Real,PES_rotor,PES_fulid}
    N::I
    B::I
    β::F
    rotor::PES_rotor
    superfulid::PES_fulid
    am2An::F
end

function C_H(ℓ,φ)
    @unpack N, B, β, rotor, superfulid, am2An = ℓ
    βE = 𝑇ᴱ_B2019(reshape(φ,3,B,N),N,B,β) - 𝑈_SuperfulidRotor(reshape(φ,3,B,N),N,B,rotor,superfulid,am2An)*β
    return βE
end

function H(ℓ,φ)
    @unpack N, B, β, rotor, superfulid, am2An = ℓ
    βE = 𝑇ᴱ_B2019(reshape(φ,3,B,N),N,B,β) - 𝑈_SuperfulidRotor(reshape(φ,3,B,N),N,B,rotor,superfulid,am2An)*β
    return βE
end

# using Pathfinder
using TransformVariables
using TransformedLogDensities: TransformedLogDensity

function (Problem::SuperfulidRotor)(φ)
    @unpack N, B, β, rotor, superfulid, am2An = Problem
    βE = 𝑇ᴱ_B2019(reshape(φ,3,B,N),N,B,β) - 
    β*𝑈_SuperfulidRotor(reshape(φ,3,B,N),N,B,rotor,superfulid,am2An)
    return βE
end

function (Problem::Classic_SuperfulidRotor)(φ)
    @unpack N, B, β, rotor, superfulid, am2An = Problem
    βE = β*𝑈_SuperfulidRotor(reshape(φ,3,1,N),N,1,rotor,superfulid,am2An)
    return -βE
end

LogDensityProblems.capabilities(::SuperfulidRotor) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.dimension(ℓ::SuperfulidRotor) = 3*N*B
LogDensityProblems.logdensity(ℓ::SuperfulidRotor, x) = H(ℓ,x)

LogDensityProblems.capabilities(::Classic_SuperfulidRotor) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.dimension(ℓ::Classic_SuperfulidRotor) = 3*N*B
LogDensityProblems.logdensity(ℓ::Classic_SuperfulidRotor, x) = C_H(ℓ,x)

# using BenchmarkTools
# LogDensityProblems.dimension(P) #Works

# @benchmark LogDensityProblems.logdensity(P, x0) #Works

# @benchmark ReverseDiff.gradient(x->LogDensityProblems.logdensity(P, x),x0)
# ∇ℓ = Hamiltonian(metric, ℓ, ReverseDiff)

dim = 3*N*B
T = as(Array, dim);
P = TransformedLogDensity(T, Problem);  
Class_P = TransformedLogDensity(as(Array, 3N), Classic_Problem)
# ∇P = ADgradient(:ForwardDiff, P);

using AdvancedMH
using Distributions
# Set up our sampler with a joint multivariate Normal proposal.
spl = MetropolisHastings(RandomWalkProposal(MvNormal(zeros(3N*B), I)))

# σ² = 0.01
# spl = MALA(x -> MvNormal((σ² / 2) .* x, σ² * I))
# model_with_ad = LogDensityProblemsAD.ADgradient(Val(:ForwardDiff), Problem)

N_n = 200
# Sample from the posterior.
chain = sample(Problem, spl, N_n;)
# E = [chain[i].lp for i in eachindex(chain)]

lP_samples = 7_000
lP_adapts = 5_000
initθ = randn(3N) .+ 2
metric = DiagEuclideanMetric(3N)
hamiltonian = Hamiltonian(metric, Class_P, ReverseDiff)
initial_ϵ = find_good_stepsize(hamiltonian, initθ)
# initial_ϵ = 0.5
integrator = Leapfrog(initial_ϵ)
proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
# proposal = HMCDA(integrator,10.0)
adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
# adaptor = MassMatrixAdaptor(metric)
# adaptor = StepSizeAdaptor(0.8, integrator)
lP, lP_stats = sample(hamiltonian, proposal, initθ, lP_samples, adaptor, lP_adapts;progress=true)

E = [-lP_stats[i][4] for i in eachindex(lP_stats)]


lP_samples = 2000
lP_adapts = 1000
initθ = begin
    tlp = lP[end]
    initp = zeros(3,B,N)
    for i in eachindex(initp[1,:,1])
        initp[:,i,:] = tlp
    end
    initp[:]
end
metric = DiagEuclideanMetric(3N*B)
hamiltonian = Hamiltonian(metric, P, ReverseDiff)
initial_ϵ = find_good_stepsize(hamiltonian, initθ)
# initial_ϵ = 0.5
integrator = Leapfrog(initial_ϵ)
proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
# proposal = HMCDA(integrator,10.0)
# adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
# adaptor = MassMatrixAdaptor(metric)
adaptor = StepSizeAdaptor(0.8, integrator)
lP, lP_stats = sample(hamiltonian, proposal, initθ, lP_samples, adaptor, lP_adapts;progress=true)

E = [-lP_stats[i][4] for i in eachindex(lP_stats)]