"""
## Rotor
`N`: Number of Particle
`B`: Beads of Path
`Z`: Atomic Number

# Set the Model
`N`: Number of Particle
`B`: Beads of Path
`T`: Temperature/K

`U`: Atomic unit (default)
`am2An` Unit Transform PE(Atomic*am2An)
"""
struct ClassicRotor{I<:Integer,F<:Real,PES_r,PES_f}
    N::I
    τ::F
    E2e::F
    rotor::PES_r
    superfluid::PES_f
end

function ClassicRotor(ℓ::SuperfluidRotor)
    @unpack N, B, τ, E2e, rotor, superfluid = ℓ
    return ClassicRotor(N,τ*B,E2e,rotor,superfluid)
end

LogDensityProblems.capabilities(::ClassicRotor) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.dimension(ℓ::ClassicRotor) = 3*N

function (Problem::ClassicRotor)(φ)
    @unpack N, τ, E2e, rotor, superfluid = Problem
    return -𝑈_SuperfluidRotor(reshape(φ,3,1,N),N,1,τ,rotor,superfluid;E2e)
end

function C2Q_init(lP,ℓ)
    @unpack N,B = ℓ
    tlp = lP[end]
    initp = zeros(3,B,N)
    for i in eachindex(initp[1,:,1])
        initp[:,i,:] = tlp
    end
    return initp[:]
end
