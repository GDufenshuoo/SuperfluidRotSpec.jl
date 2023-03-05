using MCIntegration

include("SuperfluidRotSpec.jl")
include("Bosonic.jl")

const N = 2
const B = 1
const T = 15.0
const k_b = 3.1668105084779793e-6
const β = 1/(k_b*T)
const 𝑘 = B/β
const All = N*B

X, Y, Z = Continuous(0.0, 1π, 3*All+2), Continuous(0.0, 1π, 3*All+2), Continuous(0.0, 1π, 3*All+2)
C = Dist.CompositeVar(X, Y, Z, size=3*All+2)


    ans = integrate(var=C, dof=All, neval=100_000_00, print=-1, solver=:vegasmc) do cvars, c
        x,y,z = cvars
        return -log(Wᴮ(x,y,z))/𝑘
    end

    integrate(var = Dist.CompositeVar(X, Y, Z), dof = 1, print = -1) do var, c
        X,Y,Z = var
        (X[1]^2+ Z[1]^2 + Y[1]^2 < 1.0) ? 1.0 : 0.0
    end

     integrate(var = X, dof = 3, print = -1) do X, c
        (X[1]^2 + X[2]^2 + X[3]^2 < 1.0) ? 1.0 : 0.0
    end
