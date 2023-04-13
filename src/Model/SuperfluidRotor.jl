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
    rotor::PES_r
    superfluid::PES_f
end

function SuperfluidRotor(N::Int, B::Int, RB::Int, T::Real, 
    file::String, rotor::String, superfluid::String;
    U::Unit{Float64}=Atomicᵁ, L2l::Real=1.0, E2e::Real=1.0)

    @unpack mₑ, ħ, Eᵁₖ = U
    β = 1/(Eᵁₖ*T)
    τ = β/B
    
    return SuperfluidRotor(
        N,B,RB,τ,E2e,
        set_potention(load(file)[rotor];L2l),
        set_potention(load(file)[superfluid];L2l)
        )
end

function (Problem::SuperfluidRotor)(φ)
    @unpack N, B, RB, τ, rotor, superfluid,E2e = Problem
    φ_Rotor = 5*RB
    βE = (
        𝑇ᴱ_B2019(reshape(φ[begin:φ_Rotor-1],3,B,N),N,B,τ) - 
        𝑈_SuperfluidRotor(
            reshape(φ[begin:φ_Rotor-1],3,B,N),
            reshape(φ[φ_Rotor:end],5,RB),
            N,B,RB,τ,rotor,superfluid;E2e))
    return βE
end

function 𝑈_SuperfluidRotor(x,Rx,N::Int,B::Int,RB::Int,τ::Real,rotor,superfluid;E2e=1.0)
    U1 = 0.0
    U2 = 0.0
    βU3 = β𝑈_Rotor(Rx)
    for i in 1:N
        for b in 1:B
            r = norm(x[:,b,i].-Rx[1:3,fld(b+RB-1,RB)])
            cos = ix_rot_yz(x[:,b,i],Rx[4:5,fld(b+RB-1,RB)])/r
            U1 += (r > 70.0 ? Inf : rotor(r,cos))
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

# using Plots

# const taunit = 1.4387752224
# const temprt = 0.15
# const B = 0.202_857  # Rotational constant
# nslice = 512
# τ =  taunit / (0.37 * nslice)  # Time step
# const grid_size = 100

"""
# Calculate  legendre polynomials
**Almost as fast as it can be in this structure**

Give out the Array{l} of Lₗ₋₁(x) 
"""
function legendre_polynomials(Lₘ::Int, x::Real)
    legendre = Vector{typeof(x)}(undef, Lₘ+1)
    legendre[1] = 1
    legendre[2] = x
    for i in 2:Lₘ
        ri = 1/i
        legendre[i+1] = (2-ri) * x * legendre[i] - (1-ri) * legendre[i-1]
    end
    return legendre
end

function propagator_rest(i::Int,τ::Real,B::Real)
    return exp(-τ*B*(i-1)*i)*(2i-1)/(4pi)
end

function propagator_element(Lₘ::Int, x::Real,τ::Real,B::Real)
    propagator = 0.0
    legendre = legendre_polynomials(Lₘ::Int, x::Real)
    # l = i +1 
    for i in eachindex(legendre)
        propagator += legendre[i]*propagator_rest(i,τ,B)
    end
    return abs(propagator)
end

function Rotation_x(x::Array,θ::Real)
    [
        1       0       0
        0   cos(θ) -sin(θ)
        0   sin(θ)  cos(θ)
    ]*x
end

function Rotation_y(x::Array,θ::Real)
    [
        cos(θ)  0   sin(θ)
            0   1   0
       -sin(θ)  0   cos(θ)
    ]*x
end

function Rotation_z(x::Array,θ::Real)
    [
        cos(θ) -sin(θ)  0
        sin(θ)  cos(θ)  0
        0       0       1
    ]*x
end

"""
Just Care about the x axis
"""
function ix_Rotation_y(x::Array,θ::Real)
    return (cos(θ)*x[1]+sin(θ)*x[3])
end

"""
Just Care about the x axis
"""
function ix_Rotation_z(x::Array,θ::Real)
    return (cos(θ)*x[1]-sin(θ)*x[2])
end

function ix_rot_yz(x::Array,θ::Array)
    x = ix_Rotation_y(x,θ[1])
    return x = ix_Rotation_y(x,θ[2])
end
