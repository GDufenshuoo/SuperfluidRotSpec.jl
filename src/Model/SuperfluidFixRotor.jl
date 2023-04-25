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
    rotor::PES_r
    superfluid::PES_f
end

function Change_ModelRotor(ℓ::SuperfluidFixRotor,pN::Int,pB::Int,T::Real)
    @unpack N, B, τ, rotor, superfluid = ℓ
    return SuperfluidFixRotor(pN,pB,τ/(T*B),rotor,superfluid)
end

function SuperfluidFixRotor(N::Int, B::Int, T::Real, 
    file::String, rotor::String, superfluid::String;
    U::Unit{Float64}=Atomicᵁ, L2l::Real=1.0, E2e::Real=1.0)

    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    τ = β/B
    
    return SuperfluidFixRotor(
        N,B,τ,
        set_potention(load(file)[rotor];L2l,E2e),
        set_potention(load(file)[superfluid];L2l,E2e)
        )
end

function (Problem::SuperfluidFixRotor)(φ)
    @unpack N, B, τ, rotor, superfluid = Problem
    τ⁻¹ = 1/τ
    x = reshape(φ,3,B,N)

    L_X = LX(Path_L(x),Path_X(x),τ⁻¹)
    
    βE = Tb(N,L_X)-U(x,Problem)
    
    return βE
end

function U(x,Problem::SuperfluidFixRotor)
    @unpack N, B, τ, rotor, superfluid = Problem
    U = 0.0
    for i in 1:N
        for b in 1:B
            r = sqrt(
                (x[1,b,i])^2 + 
                (x[2,b,i])^2 + 
                (x[3,b,i])^2
                )
            cos = x[1,b,i]/r
            U += rotor(r,cos)#(r > 70.0 ? Inf : )
    end end

    for i in 2:N
        for j in 1:i-1
            for b in 1:B
                r = sqrt(
                    (x[1,b,i].-x[1,b,j])^2 + 
                    (x[2,b,i].-x[2,b,j])^2 + 
                    (x[3,b,i].-x[3,b,j])^2
                    )
                U += superfluid(r)
    end end end
    return U*τ
end
