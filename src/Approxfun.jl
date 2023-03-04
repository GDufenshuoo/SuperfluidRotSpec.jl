using ApproxFun
x = Fun(-8..8)
L = -𝒟^2/2 + x^2/2
S = space(x)
B = Dirichlet(S)
λ, v = ApproxFun.eigs(B, L, 500,tolerance=1E-10)

d = (-1..1)^2
x,y = Fun(d)

