using HCubature
# Integrating cos(x) between 1.0 and 2.0
hcubature(x -> cos(x[1]), [1.0], [2.0])
# Integrating cos(x1)sin(x2) with domains of [1.0,2.0] for x1 and [1.1,3.0] for x2
hcubature(x -> cos(x[1]) * sin(x[2]), [0.0, 1.0], [0.0, 1.0])

using ApproxFun, Cubature

# 定义高维函数
f(x) = exp(-sum(x.^2))

# 定义积分变量范围和积分精度
xrange = [-1.0, 1.0]
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


using LinearAlgebra, SpecialFunctions, Plots, ApproxFun
x = Fun(identity,0..10)
y = Fun(identity,0..10)
f = sin(x^2+y^2)
g = cos(x)

h = f + g^2
r = roots(h)
rp = roots(h')

plot(h; label="f + g^2")
scatter!(r, h.(r); label="roots")
scatter!(rp, h.(rp); label="extrema")
