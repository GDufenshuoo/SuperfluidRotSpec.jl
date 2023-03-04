using MCIntegration

include("SuperfluidRotSpec.jl")
include("Bosonic.jl")

N = 2
B = 64
T = 0.15
k_b = 3.1668105084779793e-6
β = 1/(k_b*T)
𝑘 = B/β

X, Y, Z = Continuous(0.0, 1π, 3*N*B+N), Continuous(0.0, 1π, 3*N*B+N), Continuous(0.0, 1π, 3*N*B+N)
C = Dist.CompositeVar(X, Y, Z, size=3*N*B+1)

    return integrate(var=C, dof=3*N*B, neval=100000, print=-1, solver=:mcmc) do idx, cvars, c
        x,y,z = cvars
        Wᴮ(x,y,z)
    end

    return integrate(var=C, dof=1, neval=100000, print=-1, solver=:vegas) do cvars, c
        x,y,z = cvars
        Wᴮ(x,y,z)
    end

    return integrate(var=C, dof=1, neval=100000, print=-1, solver=:vegasmc) do cvars, c
        x,y,z = cvars
        Wᴮ(x,y,z)
    end

    return integrate(var=C, dof=1, neval=20, print=-1, solver=:vegas) do cvars, c
        x, y, z = cvars
        println(x[1],y[1],z[1])
        println(x[2],y[2],z[2])
        return 1.0 / (1.0 - cos(x[1]) * cos(y[1]) * cos(z[1])) / π^3
    end
    return integrate(var=C, dof=1, neval=100000, print=-1, solver=:vegasmc) do cvars, c
        x, y, z = cvars
        return 1.0 / (1.0 - cos(x[1]) * cos(y[1]) * cos(z[1])) / π^3
    end
