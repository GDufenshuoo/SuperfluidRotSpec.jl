B = 10
N = 3
x = rand(3,B,N)
𝑇ᴱ(x,N,B,1e7)
Zygote.gradient(x->𝑇ᴱ(x,N,B,1e10),x)