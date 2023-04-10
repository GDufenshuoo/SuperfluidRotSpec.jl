"""
## Electric Interaction of Atom
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
struct SuperfulidRotor{I<:Integer,F<:Real,PES_rotor,PES_fulid}
    N::I
    B::I
    β::F
    rotor::PES_rotor
    superfulid::PES_fulid
    am2An::F
end

function SuperfulidRotor(N::Int64, B::Int64, T::Float64, 
    file::String, rotor::String, superfulid::String;
    U::Unit{Float64}=Atomicᵁ, am2An = 5.29177210903e-1)

    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(T)
    
    return SuperfulidRotor(
        N,B,β,
        set_potention(load(file)[rotor]),
        set_potention(load(file)[superfulid]),
        am2An
        )
end

function (Problem::SuperfulidRotor)(φ)
    @unpack N, B, β, rotor, superfulid, am2An = Problem
    βE = 𝑇ᴱ_B2019(reshape(φ,3,B,N),N,B,β) - 
    β*𝑈_SuperfulidRotor(reshape(φ,3,B,N),N,B,rotor,superfulid,am2An)
    return -βE
end

function 𝑈_SuperfulidRotor(x,N::Int,B::Int,rotor,superfulid,am2An::Real)
    U1 = 0.0
    U2 = 0.0
    x .*= am2An
    for i in 1:N
        for b in 1:B
            r = sum(abs2,x[:,b,i])
            cos = x[1,b,i]/sqrt(r)
            U1 += (r > 30.0 ? Inf : rotor(r,cos))
    end end
    for i in 2:N
        for j in 1:i-1
            for b in 1:B
                U2 += superfulid(norm(x[:,b,i].-x[:,b,j]))
    end end end
    return (U1+U2)/B
end
