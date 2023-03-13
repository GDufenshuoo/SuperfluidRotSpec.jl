# SuperfluidRotSpec

## We can begin with An Atom Model: Li
```
N = 3
B = 2^10
T = 0.1
Z = 3
eꜛ = 2
"""Atom_Model(N::Int64, B::Int64, T::Float64, Z::Int, eꜛ::Int;U::Unit{Float64}=Atomicᵁ)"""
Problem = Atom_Model(N,B,T,Z,eꜛ)
"""runHMC(Problem,warmup::Int,Num::Int)"""
result = runHMC(Problem,3,100)
```

