function eH(Phy,K)
    E = zeros(K)
    move = similar(E)
    e = 0.0
    norm = 1.0
    Phy,E[1] = H(Phy)
    for i in 2:K
        norm /= i
        Phy,move[i],e = H(Phy)
        E[i] = E[i-1]+e*norm
    end
    # println(E)
    return Phy,move,sum(E)
end

function H(Phy)
    E = 0.0
    move = randn()
    U(x,move) = -1/abs(x) + move^2
    if exp(-(U(Phy+move,move)-U(Phy,move))) > rand()
        Phy += move
    else
        move = 0.0
    end
    return Phy,move,U(Phy,move)
end

function Solver(N,pt,K)
    P = rand(N)

end

a = eH(1.0,20)
for i in 1:100
    a = eH(a[1],10)
end
a