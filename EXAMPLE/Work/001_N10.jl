# using Pkg
# cd("SuperfluidRotSpec.jl")
# Pkg.activate(".")
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

function Savefile(file,b::String)
    open(file, "a") do io
        write(io,"\n"*">"^20*b*"\n\n")
    end
    println(">"^10*"\n Savefile $(file) \n")
end

Num = "001_single"
using Base.Threads

Wlist = [[N,B,T] for N in [10],B in [64],T in [0.37]]

const ℓ = SuperfluidFixRotor(
    1, 1, 1.0, "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2";
    L2l = 1/5.29177210903e-1,
    E2e = 1/315775.13
);#3.1668105084779793e-6/315775.13

for i in eachindex(Wlist)

    N,B,T = Wlist[i]
    N = Int(N)
    B = Int(B)

    OCSpH2 = Change_ModelRotor(ℓ,N,B,T)

    m = (15.9949146221+12.0000000+31.97207069)/5.4858e-4
    R_OCSpH2 = SuperfluidRotor(
        OCSpH2, T, m, 2, 0.202_857
    );#3.1668105084779793e-6/315775.13
    C_OCSpH2 = ClassicRotor(OCSpH2);

    println("$(Num) N$(N) B$(B) $(Int(fld(100T,1))) Class")

    ClP_w,ClP,ClPststs = runHMC(C_OCSpH2;
    lP_samples = 300, lP_adapts = 200, initθ = rand(3N), showpro=true);

    Savefile("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Class.out",[ClP_w,ClP,ClPststs])
    jldsave("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Class";ClP_w,ClP,ClPststs)
    
    println("$(Num) N$(N) B$(B) $(Int(fld(100T,1))) fix rotor")
    
    FlP_w,FlP,FlPststs = runHMC(OCSpH2;
    lP_samples = 100, lP_adapts = 50, initθ = C2Q_init(ClP,OCSpH2), showpro=true)
    Savefile("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Fix.out",[FlP_w,FlP,FlPststs])
    jldsave("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_Fix";FlP_w,FlP,FlPststs)
    Savefile("$(Num)_Output","Finished $(Num) N$(N) B$(B) $(Int(fld(100T,1))) fix rotor")

    println("$(Num) N$(N) B$(B) $(Int(fld(100T,1))) nofix rotor")

    lP_w,lP,lPststs = runHMC(R_OCSpH2;
    lP_samples = 100, lP_adapts = 50, initθ = append!(FlP[end],C2Q_init(R_OCSpH2)), showpro=true)
    Savefile("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_noFix.out",[lP_w,lP,lPststs])
    jldsave("$(Num)/file_$(Num)_N$(N)_B$(B)_$(Int(fld(100T,1)))_noFix";lP_w,lP,lPststs) 
    Savefile("$(Num)_Output","Finished $(Num) N$(N) B$(B) $(Int(fld(100T,1))) nofix rotor")

end

