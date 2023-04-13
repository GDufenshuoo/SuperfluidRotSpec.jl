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
struct SuperfluidFixRotor{I<:Integer,F<:Real,PES_r,PES_f}
    N::I
    B::I
    τ::F
    E2e::F
    rotor::PES_r
    superfluid::PES_f
end

function SuperfluidFixRotor(N::Int, B::Int, T::Real, 
    file::String, rotor::String, superfluid::String;
    U::Unit{Float64}=Atomicᵁ, L2l::Real=1.0, E2e::Real=1.0)

    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    τ = β/B
    
    return SuperfluidFixRotor(
        N,B,τ,E2e,
        set_potention(load(file)[rotor];L2l),
        set_potention(load(file)[superfluid];L2l)
        )
end

function (Problem::SuperfluidFixRotor)(φ)
    @unpack N, B, τ, rotor, superfluid,E2e = Problem
    βE = (
        𝑇ᴱ_B2019(reshape(φ,3,B,N),N,B,τ) - 
        𝑈_SuperfluidFixRotor(reshape(φ,3,B,N),N,B,τ,rotor,superfluid;E2e))
    return βE
end

function 𝑈_SuperfluidFixRotor(x,N::Int,B::Int,τ::Real,rotor,superfluid;E2e=1.0)
    U1 = 0.0
    U2 = 0.0
    for i in 1:N
        for b in 1:B
            r = norm(x[:,b,i])
            cos = x[1,b,i]/r
            U1 += (r > 70.0 ? Inf : rotor(r,cos))
    end end
    for i in 2:N
        for j in 1:i-1
            for b in 1:B
                U2 += superfluid(norm(x[:,b,i].-x[:,b,j]))
    end end end
    return (U1+U2)*τ*E2e
end
