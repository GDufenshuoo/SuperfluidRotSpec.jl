using MCIntegration

include("SuperfluidRotSpec.jl")
include("Bosonic.jl")

const N = 1
const B = 100
const T = 0.1
const k_b = 3.1668105084779793e-6
const β = 1/(k_b*T)
const 𝑘 = B/β
const A = N*B

X, Y, Z = Continuous(0.0, 1π, 3*A+2), Continuous(0.0, 1π, 3*A+2), Continuous(0.0, 1π, 3*A+2)
C = Dist.CompositeVar(X, Y, Z, size=3*A+2)

ans = integrate(var=C, dof=A, neval=1000_000, print=-1, solver=:mcmc) do idx,cvars, c
    x,y,z = cvars
    return -log(Wᴮ(x,y,z))*U(x,y,z)/β
    return U(x,y,z)
end



