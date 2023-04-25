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
    Linear_rotor::Function
    rotor::PES_r
    superfluid::PES_f
end

function SuperfluidRotor(ℓ::SuperfluidFixRotor, T, m::Real, RB::Int,LB::Real)
    @unpack N, B, τ, rotor, superfluid = ℓ
    
    return SuperfluidRotor(
        N,B,m,RB,τ,
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
        N,B,m,rRB,τ,
        LinearRotor(1.4387752224LB/(T*RB)),
        set_potention(load(file)[rotor];L2l),
        set_potention(load(file)[superfluid];L2l)
        )
end

function (Problem::SuperfluidRotor)(φ)
    @unpack N, B, m, rRB, τ, Linear_rotor, rotor, superfluid = Problem
    τ⁻¹ = 1/τ
    φm = 3N*B
    RB = fld(B,rRB)
    Rxθ = reshape((@view φ[φm+1:end]),5,RB)
    Rθ = @view Rxθ[4:5,:]
    Rx = @view Rxθ[1:3,:]

    x = reshape(φ[begin:φm],3,B,N)

    L_X = LX(Path_L(x),Path_X(x),τ⁻¹)

    βE = (
        Path(Rx)*τ*m +
        Tb(N,L_X) + 
        -U(x,Rx,Rθ,Problem)
        )
    return βE
end

function U(x,Rx,Rθ,Problem::SuperfluidRotor)
    @unpack N, B, m, rRB, τ, Linear_rotor, rotor, superfluid = Problem

    U = 0.0
    βU3 = 0.0
    rx = zeros(Real,3)

    βU3 += sum(@. log(Linear_rotor(Rθ[:,1]-Rθ[:,fld(B,rRB)])))
    for rb in 2:fld(B,rRB)
        βU3 += sum(@. log(Linear_rotor(Rθ[:,rb]-Rθ[:,rb-1])))
    end

    for i in 1:N
        for b in 1:B
            rb = fld(b+rRB-1,rRB)
            rx = (x[:,b,i].-Rx[:,rb])
            # rx = [x[1,b,i].-Rx[1,rb],x[2,b,i].-Rx[2,rb],x[3,b,i].-Rx[3,rb]]
            r = sqrt(
                (x[1,b,i])^2 + 
                (x[2,b,i])^2 + 
                (x[3,b,i])^2
                )
            cosθ = ix_rot_yz(rx,Rθ[1,rb],Rθ[2,rb])/r
            U += (r > 70.0 ? Inf : rotor(r,cosθ))
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

    return U*τ + βU3
end

