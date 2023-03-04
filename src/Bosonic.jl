"""
    Bosonic Sampling 
"""
function Wᴮ(X,Y,Z)
    E_n = Eₙ(X,Y,Z)
    W = zeros(N+1)
    W[1] = 1.0
    for n in 1:N
        for k in 1:n
            W[n+1] += exp(-E_n[k])*W[n-k+1]
        end
        W[n+1] /= n
    end
    return copy(W[N])
end

"""
    Calculate 𝑝 part
"""
function Eₙ(X,Y,Z)
    L = zeros(N)
    E_n = similar(L)
    E_x = similar(L)

    for n in 1:N, b in 2:B
        L[n] += r2(X,Y,Z,O(n,b-1),O(n,b))
    end
    for k in 2:N
        E_x[k] = r2(X,Y,Z,O(k-1,1),O(k,2))
    end
    for k in 1:N
        E_n[k] = (
            sum(E_x[(N-k+2):N]) +
            r2(X,Y,Z,O(N,1),O(N-k+1,2)) +
            sum(L[(N-k+1):N])
            )*𝑘
    end
    return E_n
end

function O(n,b)
    return (n-1)*B + b
end

function r2(X,i,j)
    return sum(abs2,X[i]-X[j])
end

function r2(X,Y,Z,i,j)
    return r2(X,i,j)+r2(Y,i,j)+r2(Z,i,j)
end
