using SuperfluidRotSpec
using JLD2
using DelimitedFiles

function Savefile(file,a)
    open(file, "w") do io
        for i in eachindex(a)
            writedlm(io,a[i]) 
            write(io,"\n"*">"^20*"\n\n")
        end
    end
    println(">"^10*"\n Savefile $(file) \n")
end

Num = "001_full"
using Base.Threads

Wlist = [[N,B,T] for N in 1:12,B in [24,64,128,256,512],T in [0.37,0.15]]

Threads.@threads for i in eachindex(Wlist)

    N,B,T = Wlist[i]
    N = Int(N)
    B = Int(B)

    OCSpH2 = SuperfluidFixRotor(
        N, B, T, "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
        L2l = 1/5.29177210903e-1,
        E2e = 1/315775.13
    );#3.1668105084779793e-6/315775.13

    m = (15.9949146221+12.0000000+31.97207069)/5.4858e-4
    R_OCSpH2 = SuperfluidRotor(
        OCSpH2, T, m, 2, 0.202_857
    );#3.1668105084779793e-6/315775.13
    C_OCSpH2 = ClassicRotor(OCSpH2);

    println("$(Num) N$(N) B$(B) $(Int(fld(100T,1))) Class")

    ClP_w,ClP,ClPststs = runHMC(C_OCSpH2;
    lP_samples = 1_000, lP_adapts = 500, initθ = rand(3N));

    Savefile("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Class.out",[ClP_w,ClP,ClPststs])
    jldsave("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Class";ClP_w,ClP,ClPststs)
    
    println("$(Num) N$(N) B$(B) $(Int(fld(100T,1))) fix rotor")

    FlP_w,FlP,FlPststs = runHMC(OCSpH2;
    lP_samples = 1_000, lP_adapts = 500, initθ = C2Q_init(ClP,OCSpH2))
    Savefile("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Fix.out",[FlP_w,FlP,FlPststs])
    jldsave("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Fix";FlP_w,FlP,FlPststs)

    println("$(Num) N$(N) B$(B) $(Int(fld(100T,1))) nofix rotor")

    lP_w,lP,lPststs = runHMC(R_OCSpH2;
    lP_samples = 1_000, lP_adapts = 500, initθ = C2Q_init(FlP,R_OCSpH2;mode=false))
    Savefile("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_noFix.out",[lP_w,lP,lPststs])
    jldsave("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_noFix";lP_w,lP,lPststs) 
end
