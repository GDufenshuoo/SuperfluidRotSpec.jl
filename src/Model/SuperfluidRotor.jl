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
    m::F
    rRB::I
    τ::F
    E2e::F
    Linear_rotor::Function
    rotor::PES_r
    superfluid::PES_f
end

function SuperfluidRotor(ℓ::SuperfluidFixRotor, T, m::Real, RB::Int,LB::Real)
    @unpack N, B, τ, E2e, rotor, superfluid = ℓ
    
    return SuperfluidRotor(
        N,B,m,RB,τ,E2e,
        LinearRotor(1.4387752224LB/(T*RB)),
        rotor,superfluid
        )
end

function SuperfluidRotor(
    N::Int, B::Int, 
    m::Real, rRB::Int,
    T::Real,LB::Real,
    file::String, rotor::String, superfluid::String;
    U::Unit{Float64}=Atomicᵁ, L2l::Real=1.0, E2e::Real=1.0)

    RB = fld(B,rRB)

    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    τ = β/B

    println("$(1.4387752224LB/(T*RB))")
    
    return SuperfluidRotor(
        N,B,m,rRB,τ,E2e,
        LinearRotor(1.4387752224LB/(T*RB)),
        set_potention(load(file)[rotor];L2l),
        set_potention(load(file)[superfluid];L2l)
        )
end

function (Problem::SuperfluidRotor)(φ)
    @unpack N, B, m, rRB, τ, Linear_rotor, rotor, superfluid,E2e = Problem
    φm = 3N*B
    RB = fld(B,rRB)
    Rxθ = reshape(φ[φm+1:end],5,RB)
    Rθ = Rxθ[4:5,:]
    Rx = Rxθ[1:3,:]
    βE = (
        β𝑇ₙ(Rx,m,RB,τ) +
        𝑇ᴱ_B2019(reshape(φ[begin:φm],3,B,N),N,B,τ) - 
        𝑈_SuperfluidRotor(
            reshape(φ[begin:φm],3,B,N),
            Rx,Rθ,
            N,B,rRB,τ,
            Linear_rotor,rotor,superfluid;E2e))
    return βE
end

function 𝑈_SuperfluidRotor(x,Rx,Rθ,N::Int,B::Int,rRB::Int,τ::Real,Linear_rotor,rotor,superfluid;E2e=1.0)
    U1 = 0.0
    U2 = 0.0
    βU3 = 0.0

    βU3 += sum(@. log(Linear_rotor(Rθ[:,1]-Rθ[:,fld(B,rRB)])))
    for rb in 2:fld(B,rRB)
        βU3 += sum(@. log(Linear_rotor(Rθ[:,rb]-Rθ[:,rb-1])))
    end

    for i in 1:N
        for b in 1:B
            rb = fld(b+rRB-1,rRB)
            rx = (x[:,b,i].-Rx[:,rb])
            r = norm(rx)
            cosθ = ix_rot_yz(rx,Rθ[:,rb])/r
            U1 += (r > 70.0 ? Inf : rotor(r,cosθ))
    end end

    for i in 2:N
        for j in 1:i-1
            for b in 1:B
                U2 += superfluid(norm(x[:,b,i].-x[:,b,j]))
    end end end

    return (U1+U2)*τ*E2e + βU3
end

function β𝑈_Rotor(θ::Array)
    U = 0.0
    for i in eachindex(θ)
        U += A(θ[i])
    end
    return U
end
