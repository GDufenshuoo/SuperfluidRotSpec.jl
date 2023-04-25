using SuperfluidRotSpec
using BenchmarkTools
using ReverseDiff
using ForwardDiff
Threads.nthreads()

N = 10
B = 512
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

P = C2Q_init([init_conf[:]],OCSpH2)
@benchmark OCSpH2(P)
@benchmark ReverseDiff.gradient(OCSpH2,P)
# @benchmark ForwardDiff.gradient(OCSpH2,P)
