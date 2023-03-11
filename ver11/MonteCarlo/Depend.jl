using UnPack
using LogDensityProblems
using Statistics
using Distributions
using LogDensityProblemsAD
using TransformedLogDensities
using TransformVariables
using DynamicHMC
using MCMCDiagnostics
import Zygote
using Random

struct ℓ
    ℓ::Function
    ∇ℓ::Function
end
