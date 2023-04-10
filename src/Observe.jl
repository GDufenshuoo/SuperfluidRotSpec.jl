<<<<<<< Updated upstream
function measure(vars, obs, weights, config) 
    # obs: prototype of the observables for each integral
        x, bin = vars #unpack the variables
        obs[1][bin[1]] += weights[1] # circle
        obs[2][bin[1]] += weights[2] # sphere
end
=======
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


function Observe(lP,l)
    @unpack N, B, β, rotor, superfulid, am2An = l
    P = size(lP,1)
    pb = zeros(2,B,N,P)
    for p in 1:P
        x = reshape(lP[p],3,B,N)*am2An
        for i in 1:N
            for b in 1:B
                r = sum(abs2,x[:,b,i])
                cos = x[1,b,i]/sqrt(r)
                pb[:,b,i,p] .= r*sqrt(1-cos^2),x[1,b,i]
    end end end
    return pb
end

# pb = Observe(lP,Problem) 
pb = Observe(lP[end-1000:end],Classic_Problem)


# using Makie
histogram2d(pb[2,1,1,:],pb[1,1,1,:],bins=(40, 20), show_empty_bins=true,
normalize=:pdf,xlim = (-5, 5),ylim = (0, 5),background_color=:Black)

histogram2d(pb[2,1,2,:],pb[1,1,2,:],bins=(40, 20), show_empty_bins=true,
normalize=:pdf,xlim = (-5, 5),ylim = (0, 5),background_color=:Black)

pythonplot()

pes = Problem.rotor

x = -5:0.005:5
y = 0.1:0.005:5
f(x, y) = begin
        x*=Classic_Problem.am2An
        y*=Classic_Problem.am2An
        r2 = x^2+y^2
        cos = x/sqrt(r2)
        E = pes[r2,cos]
        E > 1 ? 0.0 : E
    end
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
Z = map(f, X, Y)
p1 = contour(x, y, f, fill = true)
p2 = contourf(x, y, Z)
plot(p1, p2)
>>>>>>> Stashed changes
