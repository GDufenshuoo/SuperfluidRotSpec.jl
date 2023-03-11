using Random

function RandMove!(RM)
    for i in eachindex(RM)
        RM[i] = randn()
    end
end

function RandMove(io)
    RM = zeros(N)
    for i in io
        RM[i] = randn()
    end
    return RM
end

"""
# Built in Sampler
"""
function MH(𝑥,Mc,rpn)
    x = copy(𝑥)
    Acc = zeros(2)

    bin = []
    rpb = randperm(B)
    ET = Tᴬ(x).+Uᴬ(x)
    rpn = 1

    E = zeros(rpn)

    for i in 1:Mc
        io = sample(1:N, rpn, replace=false)
        for b in rpb
            # 无法组件化
            RM = RandMove!(RM)
            L = (b == 1 ? B : b-1)

            E = zeros(rpn)
            for i in eachindex(io)
                n = io[i]
                E[i] = 𝑝(x[L,n],x[b,n].+RM[i]) + 
                        -0.5*(1/abs(x[L,n]) + 1/abs(x[b,n])) +
                        -ET[b,n]
            end

            Acc[2] += 1
            if sum(x->exp(-x),E) > rand() 
                Acc[1] += 1
                for i in eachindex(io)
                    n = io[i]

                    x[b,n] += RM[i]
                    ET[b,n] += E[i]

                    push!(bin,[RM,b,io])
                end
            end
    end end
    println(Acc)
    return x,bin,ET
end

function MH(𝑥,Mc,rpn)
    x = copy(𝑥)
    Acc = zeros(2)

    bin = []
    rpb = randperm(B)
    ET = Tᴬ(x)
    EU = zeros(B,N)

    RM = zeros(N)

    for i in 1:Mc
        io = sample(1:N, rpn, replace=false)
        for b in rpb
            # 无法组件化
            L = (b == 1 ? B : b-1)
            RM = RandMove(io)
            AcR,change = δTᵒ(x,RM,b,io,ET[:,:,b])
            # AcR += -0.5*(1/abs(x[L,n]) + 1/abs(x[b,n])) -EU[b,n]
            # println(AcR)
            Acc[2] += 1
            if AcR > rand() 
                Acc[1] += 1

                x[b,:] += RM
                ET[:,:,b] += change
                push!(bin,[RM,b,io])
            end
    end end
    println(Acc)
    return x,bin,ET
end

function MH(x,Mc)
    rp = randperm(B)
    E = E(x)
    for i in 1:Mc
        io = sample(1:100, 10, replace=false)
        T = 0.0
        
    end
end

function Multilevel()
    
end