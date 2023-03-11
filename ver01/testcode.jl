include("Freemove.jl")

N = 2
B = 256

x1 = randn(B,N) .+ 3

T(pos)

Tᴬ(pos)

E(pos)

rpb = randperm(B)
T_ = Tᴬ(x)
U_ = Uᴬ(x)
rpn = 1
RM = zeros(rpn)
E = 0.0
for i in 1:Mc
    io = sample(1:N, rpn, replace=false)
    for b in rpb
        RandMove!(RM)
        E = δTᵒ(x,b,io,T)
        return RM
end end

m,n,Et = MH(x1,20000,2)

using Plots
plot(m)
histogram(abs.(m))
