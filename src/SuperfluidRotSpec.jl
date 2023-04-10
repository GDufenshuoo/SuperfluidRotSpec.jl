module SuperfluidRotSpec

<<<<<<< Updated upstream
using ApproxFun
=======
include("Units.jl")

using DynamicHMC
using LogDensityProblems
using ForwardDiff
using ReviseDiff
using UnPack


include("HMC.jl")
export runHMC

>>>>>>> Stashed changes
using FLoops
# using Tullio

greet() = print("Hello World!")

end # module SuperfluidRotSpec
