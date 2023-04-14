using SuperfluidRotSpec
using Plots
using DelimitedFiles

pythonplot()

N = 2
B = 24
T = 0.37

init_conf = zeros(3,N)

for i in eachindex(init_conf[1,:])
    init_conf[:,i] .= -10.0,8*cos(2pi*i/N),5*sin(2pi*i/N)
end

OCSpH2 = SuperfluidFixRotor(
    N, B, T, "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
    L2l = 1/5.29177210903e-1,
    E2e = 3.1668105084779793e-6);#3.1668105084779793e-6/315775.13
C_OCSpH2 = ClassicRotor(OCSpH2);
ClP_w,ClP,ClPststs = runHMC(C_OCSpH2;
lP_samples = 1_000, lP_adapts = 500, initθ = init_conf[:]);

# a_ClP = ClP

# init = [res.draws[:,i] for i in 1:size(res.draws,2)]

Op = Observe(ClP,C_OCSpH2);
histogram2d(Op[2,:,1,:][:],Op[1,:,1,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)
histogram2d(Op[2,:,:,:][:],Op[1,:,:,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)

lP_w,lP,lPststs = runHMC(OCSpH2;
lP_samples = 10, lP_adapts = 5, initθ = C2Q_init(ClP,OCSpH2))

open("output", "w") do io
    writedlm(io,[ClP_w,ClP,ClPststs])
end

Op = Observe(lP,OCSpH2)



# using Makie
histogram2d(Op[2,:,1,:][:],Op[1,:,1,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf,xlim = (-15, 15),ylim = (0, 15),background_color=:Black)

histogram2d(Op[2,:,:,:][:],Op[1,:,:,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)


# histogram2d(Op[2,1,2,:],Op[1,1,2,:],bins=(100, 100), show_empty_bins=true,
# normalize=:pdf,xlim = (-10, 10),ylim = (0, 10),background_color=:Black)

R_OCSpH2 = SuperfluidRotor(
    N, B, 12, T, 0.202_857, 
    "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
    L2l = 1/5.29177210903e-1,
    E2e = 3.1668105084779793e-6
);#3.1668105084779793e-6/315775.13

P = rand(3N*B+5*12)
R_OCSpH2(P)
P[3N*B+1:end] .= 0
lP_w,lP,lPststs = runHMC(R_OCSpH2;
lP_samples = 10, lP_adapts = 5, initθ = P)

""" 
pes = C_OCSpH2.rotor

x = -20:0.005:20
y = 0.1:0.005:20
f(x, y) = begin
        r = sqrt(x^2+y^2)
        cos = x/r
        E = pes[r,cos]
        E > 100 ? 0.0 : E
    end
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
Z = map(f, X, Y)
p1 = contour(x, y, f, fill = true)
p2 = contourf(x, y, Z)
plot(p1, p2)

"""


"""
using Pathfinder
using ReverseDiff, LogDensityProblems, LogDensityProblemsAD, TransformVariables
using TransformedLogDensities: TransformedLogDensity



function PathF(ℓ,N,B)
    println("Superfluid Rotor \n Begin to HMC")

    T = as(Array, 3*N*B)
    ℓ = TransformedLogDensity(T, ℓ)
    ∇P = ADgradient(:ReverseDiff, ℓ)
    # result_pf = pathfinder(∇P;ndraws_elbo=1000)
    multipathfinder(∇P, 10_000; nruns=10, init_scale=10)
end
"""
