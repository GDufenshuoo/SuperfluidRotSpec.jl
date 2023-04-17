using SuperfluidRotSpec
using Plots
using DelimitedFiles

pythonplot()

N = 10
B = 100
T = 0.37

OCSpH2 = SuperfluidFixRotor(
    N, B, T, "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
    L2l = 1/5.29177210903e-1,
    E2e = 1/315775.13
);#3.1668105084779793e-6/315775.13

m = (15.9949146221+12.0000000+31.97207069)/5.4858e-4
R_OCSpH2 = SuperfluidRotor(
    OCSpH2, T, m, 2, 0.202_857
);#3.1668105084779793e-6/315775.13

# m = (15.9949146221+12.0000000+31.97207069)/5.4858e-4
# R_OCSpH2 = SuperfluidRotor(
#     N, B, m, 2, T, 0.202_857, 
#     "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
#     L2l = 1/5.29177210903e-1,
#     E2e = 1/315775.13
# );#3.1668105084779793e-6/315775.13

C_OCSpH2 = ClassicRotor(OCSpH2);


ClP_w,ClP,ClPststs = runHMC(C_OCSpH2;
lP_samples = 1_000, lP_adapts = 500, initθ = rand(3*N));

FlP_w,FlP,FlPststs = runHMC(OCSpH2;
lP_samples = 1_000, lP_adapts = 500, initθ = C2Q_init(ClP,OCSpH2))

lP_w,lP,lPststs = runHMC(R_OCSpH2;
lP_samples = 1_000, lP_adapts = 500, initθ = C2Q_init(ClP,R_OCSpH2))

  p = rand(3*N)

  FOp = Observe(FlP,OCSpH2);

  lP[end][3*B*N+1:end]
  histogram(lP[end][3*B*N+1:end];bins=200)
  histogram2d(FOp[2,:,1,:][:],FOp[1,:,1,:][:],bins=(100, 100), show_empty_bins=true,
  normalize=:pdf,xlim = (-15, 15),ylim = (0, 15),background_color=:Black)
  
  histogram2d(FOp[2,:,:,:][:],FOp[1,:,:,:][:],bins=(500, 500), show_empty_bins=true,
  normalize=:pdf)

R_OCSpH2(C2Q_init([p],R_OCSpH2))

open("output_fixrotor_1", "w") do io
    writedlm(io,[FlP_w,FlP,ClPststs])
end
open("output_nofixrotor_1", "w") do io
    writedlm(io,[lP_w,lP,ClPststs])
end

R_OCSpH2(lP[1])

φm = 3N*B
RB = fld(B,2)
Rxθ = reshape(lP[1][φm+1:end],5,RB)
Rθ = Rxθ[4:5,:]
Rx = Rxθ[1:3,:]

Op = Observe(lP,R_OCSpH2);

histogram(lP[end][3*B*N+1:end];bins=200)
histogram2d(Op[2,:,1,:][:],Op[1,:,1,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf,xlim = (-15, 15),ylim = (0, 15),background_color=:Black)

histogram2d(Op[2,:,:,:][:],Op[1,:,:,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf)

f(x) = R_OCSpH2.Linear_rotor(x)
plot(-pi:0.01:pi,lf)

