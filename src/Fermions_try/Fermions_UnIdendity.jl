
"""
"""
function Samplingᵒ(px, x; Mcond=1e5)
    tx = copy(x)
    # δ = zeros(3)
    pXc = Xc(px, tx)*U(px, tx)
    for i in 1:Mcond
        for i in 1:N
            δ = randn(3).*Mᵥ/√(Mcond)
            tx[:,i] .+= δ
            X = Xc(px, tx, i)*U(px, tx) 
            if X/pXc[i] > rand()
                # println(i," ",X/pXc[i]," ",X," ",pXc[i])
                pXc[i] = X
            else
                tx[:,i] .-= δ
            end
    end end
    return tx
end

function U(px, tx)
    sum(abs2,px) + sum(abs2,tx)
end

"""Xc(px, tx)
(exp(C)*(∑exp(:::)))

123
132
213
231
312
321

123
132

123
213

123
132

"""
function Accᵒ(px, x, tx, δ, n, L; Fermions = true)
    L_C = L
    I = 0.0
    X = 0.0
    for i in 1:N
        for j in 1:N
            if i != j
                X += 𝑝(px[:,i],tx[:,j])
            else
                I += 𝑝(px[:,i],tx[:,j])
    end end end
    return exp(-𝑘*(I+L))*(1-exp(-𝑘*X/B))
end

const N = 20
const 𝑘 = 100
const Mᵥ= 1
A = rand(3)
px = zeros(3,N)
for i in 1:N
    px[:,i] = A .+ i^3*rand()
end
tx = px#rand(3,N)
Xc(px, tx)
Xc(px, tx, 1)

using BenchmarkTools

ttx = Samplingᵒ(px, tx; Mcond=1e5)
Xc(px, ttx)
Xc(px, ttx) .- Xc(px, tx)
ttx.-tx


function Xc(px, tx)
    Xc = ones(N)
    for i in 1:N
        X = 0.0
        X = 𝑝(px[:,i],tx[:,i])
        II = exp(-𝑘*X)
        for j in 1:N
            if i != j
                X = 𝑝(px[:,i],tx[:,j])
                Xc[i] *= II - exp(-𝑘*X)
    end end end
    return Xc
end

function Xc(px, tx, i)
    Xc = 1.0
    X = 1.0
    X = 𝑝(px[:,i],tx[:,i])
    II = exp(-𝑘*X)
    for j in 1:N
        if i != j
            X = 𝑝(px[:,i],tx[:,j])
            Xc *= II - exp(-𝑘*X)
    end end
    return Xc
end

"""
(exp(C)*(exp(-)-exp(X)))
(exp(O)^(B-1)*(exp(O)-exp(XX...X)))
(exp(O)^((B-1)/B)*((exp(||)-exp(X))^(∑))^(1/B))
"""
function Accᵒ(px, x, tx, δ, n, L; Fermions = true)
    L_C = L
    
    I = 0.0
    X = 0.0

    Xc = zeros(N)
    for i in 1:N
        X = 0.0
        X = 𝑝(px[:,i],tx[:,i])
        II = exp(-𝑘*X/B)
        for j in 1:N
            if i != j
                X = 𝑝(px[:,i],tx[:,j])
                Xc[i] *= II - exp(-𝑘*X/B)
            end
        end
    end



    return Xc
end

function 𝑝(A,B)
    sum(abs2,A.-B)
end
