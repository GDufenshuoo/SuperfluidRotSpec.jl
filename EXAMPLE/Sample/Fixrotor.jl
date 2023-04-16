using SuperfluidRotSpec
using Plots
using DelimitedFiles

pythonplot()

N = 3
B = 64
T = 0.37

init_conf = zeros(3,N)

for i in eachindex(init_conf[1,:])
    init_conf[:,i] .= 10.0,10*cos(2pi*i/N),10*sin(2pi*i/N)
end

OCSpH2 = SuperfluidFixRotor(
    N, B, T, "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
    L2l = 1/5.29177210903e-1,
    E2e = 1/315775.13);#3.1668105084779793e-6/315775.13
C_OCSpH2 = ClassicRotor(OCSpH2);
ClP_w,ClP,ClPststs = runHMC(C_OCSpH2;
lP_samples = 1_000, lP_adapts = 500, initθ = init_conf[:]);

lP_w,lP,lPststs = runHMC(OCSpH2;
lP_samples = 1_000, lP_adapts = 500, initθ = C2Q_init(ClP,OCSpH2));

open("output_fixrotor", "w") do io
    writedlm(io,[ClP_w,ClP,ClPststs])
end

Op = Observe(lP,OCSpH2);

histogram2d(Op[2,:,1,:][:],Op[1,:,1,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf,xlim = (-15, 15),ylim = (0, 15),background_color=:Black)

histogram2d(Op[2,:,:,:][:],Op[1,:,:,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)
