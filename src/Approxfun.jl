using ApproxFun
x = Fun(-8..8)
L = -𝒟^2/2 + x^2/2
S = space(x)
B = Dirichlet(S)
λ, v = ApproxFun.eigs(B, L, 500,tolerance=1E-10)

d = (-1..1)^2
x,y = Fun(d)

using LinearAlgebra, SpecialFunctions, Plots
using ApproxFun
x = Fun(identity,0..10)
f = sin(x^2)
g = cos(x)

using ApproxFun, Cubature

# 定义高维函数
f(x) = exp(-sum(x.^2))

# 定义积分变量范围和积分精度
xrange = (-1.0, 1.0)
npoints = 20

# 创建高维 Chebyshev 网格化对象
cheb_grid = Fourier(npoints, [xrange, xrange])

# 对高维函数进行网格化处理
f_grid = Fun(f, cheb_grid)

# 使用数值积分算法计算高维积分
result, error = hcubature(f_grid, [xrange, xrange])

# 打印积分结果和误差
println("result = ", result)
println("error = ", error)
