function Samplingᵒ(Pos_pre,Pos_now;Mcond=1e5)
    Pos_new = copy(Pos_now)
    E = Eᵒ(Pos_pre,Pos_now)
    move = randn(3).*Mᵥ

    for i in 1:Mcond
        perm = randperm(N)
        for n in 1:N
            move = randn(3).*Mᵥ
            Pos_new[:,n] .+= move
            if Accᵒ(Pos_pre, Pos_now, Pos_new, perm)#; Fermions = true)
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