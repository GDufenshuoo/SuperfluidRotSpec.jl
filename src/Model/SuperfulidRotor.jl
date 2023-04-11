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
struct SuperfluidRotor{I<:Integer,F<:Real,PES_r,PES_f}
    N::I
    B::I
    β::F
    rotor::PES_r
    superfluid::PES_f
end

function SuperfluidRotor(N::Int64, B::Int64, T::Float64, 
    file::String, rotor::String, superfluid::String;
    U::Unit{Float64}=Atomicᵁ)

    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(T)
    
    return SuperfluidRotor(
        N,B,β,
        set_potention(load(file)[rotor]),
        set_potention(load(file)[superfluid])
        )
end

function (Problem::SuperfluidRotor)(φ)
    @unpack N, B, β, rotor, superfluid = Problem
    βE = (
        𝑇ᴱ_B2019(reshape(φ,3,B,N),N,B,β) - 
        β*𝑈_SuperfluidRotor(reshape(φ,3,B,N),N,B,rotor,superfluid))
    return βE
end

function 𝑈_SuperfluidRotor(x,N::Int,B::Int,rotor,superfluid)
    U1 = 0.0
    U2 = 0.0
    for i in 1:N
        for b in 1:B
            r = norm(x[:,b,i])
            cos = x[1,b,i]/r
            U1 += (r > 30.0 ? Inf : rotor(r,cos))
    end end
    for i in 2:N
        for j in 1:i-1
            for b in 1:B
                U2 += superfluid(norm(x[:,b,i].-x[:,b,j]))
    end end end
    return (U1+U2)/B
end
