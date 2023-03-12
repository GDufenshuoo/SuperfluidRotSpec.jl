using BenchmarkTools

x = rand(3,2^10,10)
𝑈(x,10,2^10)
@benchmark 𝑈(x,10,2^10)
using ForwardDiff
using ReverceDiff
using TaylorDiff
@benchmark Zygote.gradient(x->𝑈(x,10,2^10),x)
