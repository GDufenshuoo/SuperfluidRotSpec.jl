using Random

function MH(x,Mc)
    rp = randperm(B)
    E = zeros(N,B)
    for i in 1:Mc
        io = sample(1:100, 10, replace=false)
        T = 0.0
        L = (b == 1 ? B : b-1)
        for i in io
            T += 𝑝(x[L,i],x[b,i])
        end
        return T
    end
end

function Multilevel()
    
end