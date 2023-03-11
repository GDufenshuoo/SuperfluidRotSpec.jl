
include("Depend.jl")

include("Unit.jl")
include("Model.jl")

import ForwardDiff

Problem = K_Modelly(1,1000000.0,Atomicᵁ)

dim = 3*Problem.N
T = as(Array, dim);
P = TransformedLogDensity(T, Problem);  
∇P = ADgradient(:ForwardDiff, P);

using Pathfinder

result_pf = pathfinder(∇P)

init_params = result_pf.draws[:, 1]
inv_metric = result_pf.fit_distribution.Σ
