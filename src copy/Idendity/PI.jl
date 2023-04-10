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
