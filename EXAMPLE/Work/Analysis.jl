using SuperfluidRotSpec
using JLD2
N = 10
B = 24
T = 0.37

const ℓ = SuperfluidFixRotor(
    1, 1, 1.0, "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
    L2l = 1/5.29177210903e-1,
    E2e = 1/315775.13
);#3.1668105084779793e-6/315775.13

OCSpH2 = Change_ModelRotor(ℓ,N,B,T)

m = (15.9949146221+12.0000000+31.97207069)/5.4858e-4
R_OCSpH2 = SuperfluidRotor(
    OCSpH2, T, m, 2, 0.202_857
);#3.1668105084779793e-6/315775.13
C_OCSpH2 = ClassicRotor(OCSpH2);

using Plots
Op = Observe(load("001_single/file_001_single_N10_B24_37_Fix")["FlP"],OCSpH2);
n = 10
histogram2d(Op[2,:,n,:][:],Op[1,:,n,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)
savefig("myplot_1.png") 
histogram2d(Op[2,:,:,:][:],Op[1,:,:,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)
savefig("myplot_2.png") 