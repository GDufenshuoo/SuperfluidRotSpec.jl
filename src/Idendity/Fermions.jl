"""
# The part to simulate fermions
"""
function 𝑇ᴱ(x,N,B,β)
    T = 0.0
    k = -0.5*B/β
    @floop for b in 1:B
        T += AD(x,N,B,b,k)
    end
    return -log(T^2)/β
end

"""
Under developing
"""
function AD(x,N,B,b,k)
    ans = 1.0
    A = Zygote.Buffer(zeros(),N, N)
    L = (b == 1 ? B : b-1)
    @floop for i in 1:N, j in 1:N
        A[i,j] = exp(k*𝑝(x[:,L,i],x[:,b,j]))
    end
    for j = 1:N-1
        for i = j+1:N
            if A[i,j] != 0 || A[j,j] != 0
                A[j:end,i] -= A[j,i]/A[j,j] * A[j,j:end]
    end end end
    for i in 1:N
        ans *= A[i,i]
    end
    return ans
end

