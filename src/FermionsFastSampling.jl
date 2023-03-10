
"""
这段代码是一个实现了一些采样、计算能量和接受拒绝的函数。
它们被用于某些统计物理模型的模拟，其中多粒子系统的空间位置被统计。
下面是各函数的含义和作用：

Samplingᵒ(Pos_pre,Pos_now;Mcond=1e5): 
    该函数模拟从 Pos_now 开始通过在参数空间中生成随机步长来进行 MCMC 采样。
    每个步骤从当前位置向随机方向移动。该函数使用 Eᵒ() 函数计算界面上一个二体相互作用的能量，
    同时接受或拒绝随机步骤。函数返回新的坐标 Pos_new，通过 Pos_now 向每个粒子应用随机步长计算而来。
Accᵒ(Pos_pre,Pos_now,Pos_new,perm;Fermions=true): 
    该函数计算 Metropolis-Hastings 接受拒绝准则，并且在满足条件时接受新的位向 Pos_new。调用时，需要传递先前的位置 Pos_pre，当前位置 Pos_now，新的位置 Pos_new，以及排列 perm（粒子之间的对应关系）。默认情况下，Fermions 参数是正确的，以便在处理费米子时正确地计算矩阵元素。该函数返回一个布尔值，指示是否接受新的状态。
Accᵒ(Pos_pre, Pos_now, perm): 
    该函数计算 Metropolis-Hastings 接受拒绝准则，并在满足条件时接受新的位向 Pos_new。这个函数与上面那个函数有一些区别，因为这个函数只计算一个方向的运动，并且在内部生成基于原始位置和一个随机向量变化的新的位置向量。然后使用 Accᵒ(Pos_pre,Pos_now,Pos_new,perm;Fermions=true) 中的 Pos_new 和其他参数计算接受拒绝准则。
Wᶠᵒ(Pos_pre, Pos_now, perm): 
    该函数计算在两个位置 Pos_now 和 Pos_now 之间的正向权重。这个函数计算在这两个位置上的能量差异量，并将其除以二，这是由 Metropolis-Hastings 接受拒绝准则中给出的比例常数。
Wᴮᵒ(Pos_pre, Pos_now, perm): 
    该函数计算在两个位置 Pos_now 和 Pos_now 之间的反向权重。
Eᵒ(Pos_pre, Pos_now, perm): 
    该函数计算在两个位置 Pos_now 和 Pos_now 之间的能量，这是在欧几里得空间中计算两个向量之间的距离（也就是 L2 范数）的简单计算。
Eᵒ(Pos_pre, Pos_now): 
    该函数计算在两个位置 Pos_now 和 Pos_now 之间的能量（即角标 perm 已被随机打乱或忽略）。这样计算能量，用于确定在没有重排粒子的情况下的能量。
"""
function Samplingᵒ(Pos_pre,Pos_now;Mcond=1e5)
    Pos_new = copy(Pos_now)
    E = Eᵒ(Pos_pre,Pos_now)
    move = randn(3).*Mᵥ

    for i in 1:Mcond
        perm = 
        for n in 1:N
            move = randn(3).*Mᵥ
            Pos_new[:,n] .+= move
            if Accᵒ(Pos_pre, Pos_now, Pos_new,perm)#; Fermions = true)
                print("O",move)
            else
                Pos_new[:,n] .-= move
    end end end
    return Pos_new
end

𝑘 = 1
N = 2

Pos_pre = 10*rand(3,2)
Pos_now = 10*rand(3,2)

Mᵥ = 1
a = Samplingᵒ(Pos_pre,Pos_now;Mcond=1e3).-Pos_now

Check_On = true

if Check_On
    for i in 1:100
        Pos_pre = rand(3,2)
        Pos_now = rand(3,2)
        perm = randperm(2)
        Wᴮᵒ(Pos_pre, Pos_now, perm)
        exp(-Eᵒ(Pos_pre, Pos_now))
        Wᶠᵒ(Pos_pre, Pos_now, perm)/exp(-Eᵒ(Pos_pre, Pos_now))
        println(Accᶠᵒ(Pos_pre, Pos_now, perm))
        # println(Accᴮᵒ(Pos_pre, Pos_now, perm))
    end
end

function Accᵒ(Pos_pre, Pos_now, Pos_new, perm; Fermions = true)
    println((
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new))) - 
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new, perm)))
        )/2
        )
    return (
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new))) - 
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new, perm)))
        )/2/(exp(-𝑘*Eᵒ(Pos_pre, Pos_now))
        ) > rand()
end

function Accᵒ(Pos_pre, Pos_now, perm)
    println((
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new))) - 
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new, perm)))
        )/2
        )
    return (
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new))) - 
        exp(-𝑘*(Eᵒ(Pos_pre, Pos_new, perm)))
        )/2/(exp(-𝑘*Eᵒ(Pos_pre, Pos_now))
        ) > rand()
end

function Wᶠᵒ(Pos_pre, Pos_now, perm)
    return (
        exp(-𝑘*Eᵒ(Pos_pre, Pos_now)) - 
        exp(-𝑘*Eᵒ(Pos_pre, Pos_now, perm))
        )/2
end

function Wᴮᵒ(Pos_pre, Pos_now, perm)
    return (
        exp(-𝑘*Eᵒ(Pos_pre, Pos_now)) + 
        exp(-𝑘*Eᵒ(Pos_pre, Pos_now, perm))
        )/2
end

function Eᵒ(Pos_pre, Pos_now, perm)
    E = 0.0
    for i in 1:N
        E += sum(Pos_pre[:,i]-Pos_now[:,perm[i]])
    end
    return E
end

function Eᵒ(Pos_pre, Pos_now)
    E = 0.0
    for i in 1:N
        E += sum(Pos_pre[:,i]-Pos_now[:,i])
    end
    return E
end

