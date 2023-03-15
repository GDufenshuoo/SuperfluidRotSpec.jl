using SuperfluidRotSpec

N = 1
B = 2^7
T = 0.01
Z = 1
eꜛ = 1
"""Atom_Model(N::Int64, B::Int64, T::Float64, Z::Int, eꜛ::Int;U::Unit{Float64}=Atomicᵁ)"""
Problem = Atom_Model(N,B,Z,eꜛ,T)
using BenchmarkTools
Problem(rand(3*N*B))
# @benchmark Problem(rand(3*N*B))
result = runHMC(Problem,6,1000)
# @benchmark Zygote.gradient(x->Problem(x),rand(3*N*B))
# println(show(Problem),"\n Begin to HMC")
# ℓ_dims = 3*N*B
# T = as(Array, ℓ_dims);
# ℓ = TransformedLogDensity(T, Problem);  
# ∇ℓ = ADgradient(:Zygote, ℓ);
# rng = Random.GLOBAL_RNG

# mcmc_keep_warmup(rng, ∇ℓ, 100;)
