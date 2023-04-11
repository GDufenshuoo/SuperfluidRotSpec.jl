using SuperfluidRotSpec

N = 2
B = 128
T = 0.15

OCSpH2 = SuperfluidRotor(N, B, T, "PES/OCS.PES", "OCS_paraH2", "paraH2_paraH2")
C_OCSpH2 = ClassicRotor(OCSpH2)
ClP_w,ClP,ClPststs = runHMC(C_OCSpH2 ;
lP_samples = 1_000, lP_adapts = 500, initθ = randn(3N) .+ 2)
lP_w,lP,lPststs = runHMC(OCSpH2 ;
lP_samples = 1000, lP_adapts = 500, initθ = C2Q_init(lP,OCSpH2))

Op = Observe(lP,OCSpH2)

using Plots
pythonplot()

# using Makie
histogram2d(Op[2,:,1,:][:],Op[1,:,1,:][:],bins=(100, 100), show_empty_bins=true,
normalize=:pdf,xlim = (-5, 5),ylim = (0, 5),background_color=:Black)



histogram2d(Op[2,1,2,:],Op[1,1,2,:],bins=(100, 20), show_empty_bins=true,
normalize=:pdf,xlim = (-5, 5),ylim = (0, 5),background_color=:Black)



"""pes = C_OCSpH2.rotor

x = -5:0.005:5
y = 0.1:0.005:5
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
plot(p1, p2)"""

