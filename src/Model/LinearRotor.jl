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

function Rotation_x(x,θ::Real)
    [
        1       0       0
        0   cos(θ) -sin(θ)
        0   sin(θ)  cos(θ)
    ]*x
end

function Rotation_y(x,θ::Real)
    [
        cos(θ)  0   sin(θ)
            0   1   0
       -sin(θ)  0   cos(θ)
    ]*x
end

function Rotation_z(x,θ::Real)
    [
        cos(θ) -sin(θ)  0
        sin(θ)  cos(θ)  0
        0       0       1
    ]*x
end

"""
Just Care about the x axis
"""
function ix_Rotation_y(x,θ::Real)
    return (cos(θ)*x[1]+sin(θ)*x[3])
end

"""
Just Care about the x axis
"""
function ix_Rotation_z(x,θ::Real)
    return (cos(θ)*x[1]-sin(θ)*x[2])
end

function ix_rot_yz(x,θ,ϕ)
    Ry_x = Rotation_y(x,θ)
    return ix_Rotation_z(Ry_x,ϕ)
end

"""
0.003<τB<0.6  maxerr<0.01
"""
function LinearRotor(τB)
    h = (((0.07957649676528922 / τB) + (τB * 0.005536640332578352)) + 0.026534110958300494)
    σ = (((0.2508712935565096 / τB) + -0.09420101021905938) ^ 0.4995735458437035)
    return LinearRotor(θ) = h*exp(-(σ*θ)^2)
end
