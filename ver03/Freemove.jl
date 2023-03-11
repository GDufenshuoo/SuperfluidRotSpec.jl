"""
### 用 Slater 行列式来解自由费米粒子
"""
function 𝑇ᴱ(x,N,B,β)
    T = 0.0
    k = B/β
    for b in 1:B
        T += det(AD(x,N,B,b,k))
    end
    println(T)
    return -log(T^2)/β
end

function AD(x,N,B,b,k)
    A = Zygote.Buffer(zeros(),N, N)
    L = (b == 1 ? B : b-1)
    for i in 1:N, j in 1:N
        A[i,j] = exp(-0.5k*
        𝑝(x[:,L,i],x[:,b,j]) - 
        0.5/k*𝑈(x[:,L,i],x[:,b,j])
        )
    end
    return copy(A)
end

"""
### get the 𝑝
"""
function 𝑝(a,b)
    p = 0.0
    for i in 1:3
        p += (a[i]-b[i])^2
    end
    return p
end

function 𝑈(a,b)
    U = 0.0
    U += 3/(norm(a)+1e-10)
    U += 3/(norm(b)+1e-10)
    U -= 1/(norm(a.-b)+1e-10)
    return U/B
end

function 𝑈(x,N,B)
    U = 0.0
    for i in 1:N
        for b in 1:B
            U += 3/(norm(x[:,b,i])+1e-10)
    end end
    for i in 2:N
        for j in 1:i, b in 1:B
            U -= 1/(norm(x[:,b,i].-x[:,b,j])+1e-10)
    end end
    return U/B
end