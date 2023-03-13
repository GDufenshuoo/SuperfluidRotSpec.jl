"""
# The part to simulate fermions
"""
function 𝑇ᴱ(x,N,B,β)
    T = 0.0
    k = B/β
    for b in 1:B
        T += det(AD(x,N,B,b,k))
    end
    return -log(T^2)/β
end

"""
Under developing
"""
function AD(x,N,B,b,k)
    A = Zygote.Buffer(zeros(),N, N)
    L = (b == 1 ? B : b-1)
    for i in 1:N, j in 1:N
        A[i,j] = exp(-0.5*k*𝑝(x[:,L,i],x[:,b,j]))
    end
    return copy(A)
end
