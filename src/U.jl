using LinearAlgebra

function O(n,b)
    return (n-1)*B + b
end

function Com(X,Y,Z,i,j)
    return norm([X[i]-X[j], Y[i]-Y[j], Z[i]-Z[j]])
end

function r2(X,i)
    return sum(abs2,X[i])
end

function r2(X,i,j)
    return sum(abs2,X[i]-X[j])
end

function r2(X,Y,Z,i)
    return r2(X,i)+r2(Y,i)+r2(Z,i)
end

function r2(X,Y,Z,i,j)
    return r2(X,i,j)+r2(Y,i,j)+r2(Z,i,j)
end

function U(X,Y,Z)
    L = 0.0
    for i in 2:N, j in 1:i, b in 1:B
        L -= Com(X,Y,Z,O(i,b),O(j,b))
    end
    L
end