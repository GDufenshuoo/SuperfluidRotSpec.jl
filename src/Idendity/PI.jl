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

function β𝑇ₙ(x,m::Real,B::Int,τ::Real)
    T = 0.0
    k = -1/(2τ)
    T += 𝑝(x[:,B],x[:,1])
    for b in 2:B
        T += 𝑝(x[:,b-1],x[:,b])
    end
    return m*k*T
end
