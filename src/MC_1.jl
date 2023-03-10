using MonteCarloMeasurements, StatsPlots, NNlib, Measurements, KernelDensity
using Measurements: value, uncertainty
default(lab="")
N    = 20 # Number of particles
f    = x -> σ(12x-6) # Nonlinear function
l    = 0 # Left boundary
r    = 1 # Right boundary
d    = Normal(0.5, 0.15) # The probability density of the input
m    = Measurements.:(±)(d.μ, d.σ) # For comparison to Measurements.jl using linear uncertainty propagation
my   = h(m) # output measurement
dm   = Normal(value(my), uncertainty(my)) # Output density according to Measurements
x    = Particles(N,d, permute=false) # Create particles distributed according to d, sort for visualization
y    = h(x).particles # corresponding output particles
x    = x.particles # extract vector to plot manually
xr   = LinRange(l,r,100) # x values for plotting
noll = zeros(N)


plot(f, l, r, legend=:right, xlims=(l,r), ylims=(l,r), axis=false, grid=false, lab="h(x)", xlabel="Input space", ylabel="Output space")
plot!(x->0.2pdf(d,x),l,r, lab="Input dens.")
# Estimate the true output density using a large sample
kdt = kde(f.(rand(d,100000)), npoints=200, bandwidth=0.08)
plot!(l .+ 0.2kdt.density, kdt.x, lab="True output dens.")

# This is the output density as approximated by linear uncertainty propagation
plot!(l .+ 0.2pdf.(Ref(dm),xr), xr, lab="Linear Gaussian propagation")

# Estimate the output density corresponding to the particles
kd = kde(y, npoints=200, bandwidth=0.08)
plot!(l .+ 0.2kd.density, kd.x, lab="Particle kernel dens. est.", l=:dash)

# Draw helper lines that show how particles are transformed from input space to output space
plot!([x x][1:2:end,:]', [noll y][1:2:end,:]', l=(:black, :arrow, :dash, 0.1))
plot!([x fill(l,N).+0.02][1:2:end,:]', [y y][1:2:end,:]', l=(:black, :arrow, :dash, 0.1))

# Plot the particles
scatter!(x, 0y, lab="Input particles")
scatter!(fill(l,N) .+ 0.02, y, lab="Output particles")

# Draw mean lines, these show hoe the mean is transformed using linear uncertainty propagation
plot!([d.μ,d.μ], [0,h(d.μ)], l=(:red, :dash, 0.2))
plot!([l,d.μ], [h(d.μ),h(d.μ)], l=(:red, :dash, 0.2))



using Pkg
Pkg.add.(
      ["Distributions"
    , "MonteCarloMeasurements"
    , "StatsFuns"
    , "LaTeXStrings"
    , "Plots"])


using Distributions
using MonteCarloMeasurements
using StatsFuns

# Just some setup so inequalities propagate through particles
for rel in [<,>,<=,>=]
    register_primitive(rel)
end


# function fromObs(x,y) 
    function logp(α,β)
        ℓ = 0.0
        ℓ += logpdf(Normal(0,1), α)
        ℓ += logpdf(Normal(0,2), β)
        yhat = α .+ β .* x
        ℓ += sum(logpdf.(Normal.(yhat, 1), y) )
        ℓ
    end
# end


drawcat(ℓ, k) = [argmax(ℓ + Particles(1000,Gumbel())) for j in 1:k]

asmatrix(ps...) = Matrix([ps...])'

# Kish's effective sample size,
# $n _ { \mathrm { eff } } = \frac { \left( \sum _ { i = 1 } ^ { n } w _ { i } \right) ^ { 2 } } { \sum _ { i = 1 } ^ { n } w _ { i } ^ { 2 } }$

function n_eff(ℓ)
    logw = ℓ.particles
    exp(2 * logsumexp(logw) - logsumexp(2 .* logw))
end


function h(a,b)
    # generate data
    x = rand(Normal(),100)
    yhat = a .+ b .* x
    y = rand.(Normal.(yhat, 1))
    # generate p
    logp = fromObs(x,y)

    runInference(x,y,logp)
end

# function runInference(x,y,logp)
    x = 3
    y = 4
    N = 1000 

    # initialize q
    q = MvNormal(2,100000.0) # Really this would be fit from a sample from the prior
    α,β = Particles(N,q)
    m = asmatrix(α,β)
    ℓ = sum(logp(α,β)) - Particles(logpdf(q,m))

    numiters = 60
    elbo = Vector{Float64}(undef, numiters)
    for j in 1:numiters
        α,β = Particles(N,q)
        m = asmatrix(α,β)
        ℓ = logp(α,β) - Particles(logpdf(q,m))
        elbo[j] = mean(ℓ) 
        ss = suffstats(MvNormal, m,  exp(ℓ - maximum(ℓ)).particles .+ 1/N)
        q = fit_mle(MvNormal, ss)
    end
    (α,β,q,ℓ,elbo)
# end


(α,β,q,ℓ,elbo) = h(3,4)

using LaTeXStrings
using Plots
plot(1:60, -elbo
    , xlabel="Iteration"
    , ylabel="Negative ELBO"
    , legend=false
    , yscale=:log10)
xticks!([0,20,40,60], [L"0",L"20", L"40",L"60"])
yticks!(10 .^ [3,6,9,12], [L"10^3", L"10^6",L"10^9",L"10^{12}"])
savefig("neg-elbo.svg")






