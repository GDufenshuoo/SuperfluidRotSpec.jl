using SuperfluidRotSpec
using JLD2
N = 10
B = 128
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
    Op = Observe(load("EXAMPLE/Work/file_N10_B256_37_nofix")["lP"],R_OCSpH2);
histogram2d(Op[2,:,:,:][:],Op[1,:,:,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)
savefig("myplot_2_$(N)_$(B).png") 
for n in 1:N
    histogram2d(Op[2,:,n,:][:],Op[1,:,n,:][:],bins=(100, 100), show_empty_bins=true,
    normalize=:pdf)
    savefig("myplot_1_$(N)_$(B)_$n.png") 
end
