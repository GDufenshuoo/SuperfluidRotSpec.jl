"""
# Bosonic 𝑇 using The method from the 
    DOI: 10.1073/pnas.1913365116
    "Path integral molecular dynamics for bosons"
### Almost as fast as it can be
"""
function 𝑇ᴱ_B2019(x,N,B,β)
    𝑘 = -0.5*B/β
    A = Zygote.Buffer(zeros(),N,3)

    @floop for i in 1:N
        for b in 2:B
            A[i,3] += 𝑝(x[:,b-1,i],x[:,b,1])
    end end

    # Can be better but i don't want to :3
    @floop for k in 1:N
        for i in N-k+1:N
        R = (i == N) ? N-k+1 : i
        A[i,2] += 𝑝(x[:,B,i],x[:,1,R]) + 
                    + A[i,3]
    end end

    # EN(k) = A[i,2]
    # exp(k*VB(N-k)) = A[i-k,1]
    for i in 1:N
        for k in 1:i
            if i == k
                A[i,1] += (exp(𝑘*(
                    A[i,2]
                )))/i
            else
                A[i,1] += A[i-k,1]*(exp(𝑘*(
                    A[i,2]
                )))/i
    end end end

    return -log.(copy(A[N,1]))
end



