# 

"""
    Bosonic Sampling 
"""
function Wᴮ(P)
    E_n = Eₙ(P)
    W = Zygote.Buffer(zeros(N+1),N+1)
    W[1] = 1.0
    for n in 1:N
        @floop for k in 1:n
            W[n+1] += exp(-E_n[k])*W[n-k+1]
        end
        W[n+1] /= n
    end
    return abs(copy(W[N]))
end

function Eₙ(P)
    E_n = Zygote.Buffer(zeros(N),N)
    E_x = Zygote.Buffer(zeros(N),N)
    L = Zygote.Buffer(zeros(N),N)

    @floop for n in 1:N, b in 2:B
        L[n] += sum(abs2,P[:,b,n].-P[:,b-1,n])
    end
    @floop for k in 2:N
        E_x[k] = sum(abs2,P[:,1,k-1].-P[:,2,k])
    end
    @floop for k in 1:N
        E_n[k] = (sum(E_x[(N-k+2):N]) +
                 sum(abs2,P[:,1,N].-P[:,2,N-k+1]) +
                 sum(L[(N-k+1):N]))*k
    end
    return E_n
end



# 

"""
    Bosonic Sampling 
"""
function WᴮO(P)
    E_n = EₙO(P)
    W = Zygote.Buffer(zeros(N+1),N+1)
    W[1] = 1.0
    for n in 1:N
        @floop for k in 1:n
            W[n+1] += exp(-E_n[k])*W[n-k+1]
        end
        W[n+1] /= n
    end
    return abs(copy(W[N]))
end

function EₙO(P)
    E_n = Zygote.Buffer(zeros(N),N)
    E_x = Zygote.Buffer(zeros(N),N)

    @floop for k in 2:N
        E_x[k] = sum(abs2,P[:,1,k-1].-P[:,2,k])
    end
    @floop for k in 1:N
        E_n[k] = (sum(E_x[(N-k+2):N]) +
                 sum(abs2,P[:,1,N].-P[:,2,N-k+1]))*k
    end
    return E_n
end

