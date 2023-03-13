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
struct SuperfulidRotor{I<:Integer,F<:Real,PES}
    N::I
    B::I
    β::F
    rotor::PES
    superfulid::PES
end

function SuperfulidRotor(N::Int64, B::Int64, T::Float64, 
    file::String, rotor::String, superfulid::String;
    U::Unit{Float64}=Atomicᵁ, am2An = 5.29177210903e-1)

    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    
    return SuperfulidRotor(
        N,B,β,
        set_potention(load(file)[rotor]),
        set_potention(load(file)[superfulid])
            )
end

function (Problem::SuperfulidRotor)(φ)
    @unpack N, B, β = Problem
    E = 𝑇ᴱ(reshape(φ,3,B,N),Problem::SuperfulidRotor) + 
        𝑈(reshape(φ,3,B,N),Problem::SuperfulidRotor)
    return -E
end

"""
# The part to simulate fermions
"""
function 𝑇ᴱ(x,Problem::SuperfulidRotor)
    @unpack N,B,β = Problem
    return 𝑇ᴱ_B2019(x,N,B,β)
end

function 𝑈(x,Problem::SuperfulidRotor)
    @unpack N,B,rotor,superfulid = Problem
    U1 = 0.0
    U2 = 0.0
    @floop for i in 1:N
        for b in 1:B
            r = norm(x[:,b,i])*am2An
            cos = x[1,b,i]/r
            @reduce U1 += (r > 30.0 ? 0.0 : rotor(r,cos))
    end end
    @floop for i in 2:N
        for j in 1:i
            for b in 1:B
                @reduce U2 += superfulid(norm(x[:,b,i].-x[:,b,j])*am2An)
    end end end
    return (U1+U2)/B
end
