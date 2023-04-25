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
    τB::F
    rotor::PES_r
    superfluid::PES_f
end

function ClassicRotor(ℓ::SuperfluidFixRotor)
    @unpack N, B, τ, rotor, superfluid = ℓ
    return ClassicRotor(N,τ*B,rotor,superfluid)
end

LogDensityProblems.capabilities(::ClassicRotor) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.dimension(ℓ::ClassicRotor) = 3*N

function (Problem::ClassicRotor)(φ)
    @unpack N, τB, rotor, superfluid = Problem
    return -𝑈_SuperfluidFixRotor(reshape(φ,3,1,N),N,1,τB,rotor,superfluid)
end

function C2Q_init(lP,ℓ::SuperfluidFixRotor)
    @unpack N,B = ℓ
    tlp = lP[end]
    initp = zeros(3,B,N)
    for i in eachindex(initp[1,:,1])
        initp[:,i,:] = tlp
    end
    return initp[:]
end

function C2Q_init(ℓ::SuperfluidRotor)
    @unpack N,B,rRB = ℓ
    RB = fld(B,rRB)
    return 0.001*randn(5*RB)
end
