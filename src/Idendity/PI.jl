function Path(x)
    A = 0.0
    B = size(x,2)
    N = size(x,3)
    for n in 1:N
        for (i,b) in zip(n,2:B)
            A += (x[1,b-1,i]-x[1,b,i])^2
            A += (x[2,b-1,i]-x[2,b,i])^2
            A += (x[3,b-1,i]-x[3,b,i])^2
        end
    end
    return A
end