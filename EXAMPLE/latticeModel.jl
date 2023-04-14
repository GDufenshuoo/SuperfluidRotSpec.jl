# test_code
using LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics
using Parameters, Statistics, Random
using Zygote

#Both value and gradient in same calculation:
function value_and_gradient(f, x...)
    value, back = Zygote.pullback(f, x...)
    grad = back(1)[1]
    return value, grad
end

#### Define the problem ####
struct LogNormalTest
    n::Int
end

(ℓ::LogNormalTest)(x) = -sum(x.^2)

LogDensityProblems.capabilities(::LogNormalTest) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.dimension(ℓ::LogNormalTest) = ℓ.n
LogDensityProblems.logdensity(ℓ::LogNormalTest, x) = ℓ(x)
LogDensityProblems.logdensity_and_gradient(ℓ::LogNormalTest, x) = value_and_gradient(ℓ, x)

#### Testing the problem ####
n = 100 #Should, ideally, be much larger...
ℓ = LogNormalTest(n)

x0 = randn(n)
LogDensityProblems.dimension(ℓ) #Works
LogDensityProblems.logdensity(ℓ, x0) #Works
LogDensityProblems.logdensity_and_gradient(ℓ, x0) #Works

#### HMC ####
rng = Random.GLOBAL_RNG #Works with this, but only on the CPU..?
# rng = CURAND.RNG() #Fails with this... But more is probably needed?
results = mcmc_with_warmup(rng, ℓ, 10) #Errors with "KernelError: passing and using non-bitstype argument"
