x = rand(3,2^10,10)
𝑇ᴱ(x,10,2^10,1e8)
Zygote.gradient(x->𝑇ᴱ(x,10,2^10,1e8),x)
using BenchmarkTools
@benchmark 𝑇ᴱ(x,10,2^10,1e8)
@benchmark Zygote.gradient(x->𝑇ᴱ(x,10,2^10,1e8),x)


