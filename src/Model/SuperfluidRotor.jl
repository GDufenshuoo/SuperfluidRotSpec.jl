include("LinearRotor.jl")

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
    RB::I
    τ::F
    E2e::F
    Linear_rotor::var"#LinearRotor#5"{F, F}
    rotor::PES_r
    superfluid::PES_f
end

function SuperfluidRotor(ℓ::SuperfluidFixRotor, RB::Int,LB::Real)
    @unpack N, B, τ, E2e, rotor, superfluid = ℓ
    
    return SuperfluidRotor(
        N,B,RB,τ,E2e,
        LinearRotor(1.4387752224LB/(T*RB)),
        set_potention(load(file)[rotor];L2l),
        set_potention(load(file)[superfluid];L2l)
        )
end

function SuperfluidRotor(N::Int, B::Int, RB::Int,
    T::Real,LB::Real,
    file::String, rotor::String, superfluid::String;
    U::Unit{Float64}=Atomicᵁ, L2l::Real=1.0, E2e::Real=1.0)


    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    τ = β/B
    
    return SuperfluidRotor(
        N,B,RB,τ,E2e,
        LinearRotor(1.4387752224LB/(T*RB)),
        set_potention(load(file)[rotor];L2l),
        set_potention(load(file)[superfluid];L2l)
        )
end

function (Problem::SuperfluidRotor)(φ)
    @unpack N, B, RB, τ, Linear_rotor, rotor, superfluid,E2e = Problem
    φm = 3N*B
    βE = (
        𝑇ᴱ_B2019(reshape(φ[begin:φm],3,B,N),N,B,τ) - 
        𝑈_SuperfluidRotor(
            reshape(φ[begin:φm],3,B,N),
            reshape(φ[φm+1:end],5,RB),
            N,B,RB,τ,
            Linear_rotor,rotor,superfluid;E2e))
    return βE
end

function 𝑈_SuperfluidRotor(x,Rx,N::Int,B::Int,RB::Int,τ::Real,Linear_rotor,rotor,superfluid;E2e=1.0)
    U1 = 0.0
    U2 = 0.0
    expβU3 = 0.0
    for i in 1:N
        for b in 1:B
            rx = (x[:,b,i].-Rx[1:3,fld(b+RB-1,RB)])
            r = norm(rx)
            Urx = rx/r
            # println("Urx $Urx $r ")
            # println("$(ix_rot_yz(x[:,b,i],Rx[4:5,fld(b+RB-1,RB)]))")
            cosθ = ix_rot_yz(Urx,Rx[4:5,fld(b+RB-1,RB)])
            # println("rx $rx r $r")
            # println("cos $cos rx[1]/r $(rx[1]/r)")
            θ = acos(cosθ) - acos(Urx[1])
            expβU3 *= Linear_rotor(θ)
            U1 += (r > 70.0 ? Inf : rotor(r,cosθ))
    end end
    for i in 2:N
        for j in 1:i-1
            for b in 1:B
                U2 += superfluid(norm(x[:,b,i].-x[:,b,j]))
    end end end
    return (U1+U2)*τ*E2e + log(expβU3)
end

function β𝑈_Rotor(θ::Array)
    U = 0.0
    for i in eachindex(θ)
        U += A(θ[i])
    end
    return U
end
