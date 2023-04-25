
struct LX
    E_::Array
    Ex_::Matrix
    τ⁻¹::Real
end

"""
# Bosonic 𝑇 using The method from the 
    DOI: 10.1073/pnas.1913365116
    "Path integral molecular dynamics for bosons"
"""
function Tb(N,LX)
    # A = zeros(Real,BigFloat,N)
    A = zeros(Real,N)
    for i in 1:N
        ExV!(A,i,LX)
    end
    return A[end]
end

function ExV!(A,N,LX)
    # eV = BigFloat(0.0)
    eV = (0.0)
    if N == 0
        return 1.0
    end
    for k in 1:N-1
        eV += A[N-k]*ExE(k,N,LX)
    end
    eV += ExE(N,N,LX)
    A[N] = eV/N
end

function ExV(N,LX)
    # eV = BigFloat(0.0)
    eV = (0.0)
    if N == 0
        return 1.0
    end
    for k in 1:N
        eV += ExV(N-k,LX)*ExE(k,N,LX)
    end

    return eV/N
end

function ExE(k,N,LX)
    @unpack E_,Ex_,τ⁻¹ = LX
    # E = BigFloat(0.0)
    E = (0.0)
    for l in N-k+1:N-1
        E += E_[l] + Ex_[l,l+1]
    end
    E += E_[N] + Ex_[N,N-k+1]


    return exp(-τ⁻¹*E)
end

function Path_L(x)
    B = size(x,2)
    N = size(x,3)
    A = zeros(Real,N)

    for n in 1:N
        for (i,b) in zip(n,2:B)
            A[i] += (x[1,b-1,i]-x[1,b,i])^2
            A[i] += (x[2,b-1,i]-x[2,b,i])^2
            A[i] += (x[3,b-1,i]-x[3,b,i])^2
        end
    end
    return A
end

function Path_X(x)
    B = size(x,2)
    N = size(x,3)
    A = zeros(Real,N,N)

    for j in 1:N
        for i in 1:j
            A[i,j] += (x[1,1,i]-x[1,B,j])^2
            A[i,j] += (x[2,1,i]-x[2,B,j])^2
            A[i,j] += (x[3,1,i]-x[3,B,j])^2
    end end
    return A
end

function Path_L(A,x)
    B = size(x,2)
    N = size(x,3)
    for n in 1:N
        for (i,b) in zip(n,2:B)
            A[i] += (x[1,b-1,i]-x[1,b,i])^2
            A[i] += (x[2,b-1,i]-x[2,b,i])^2
            A[i] += (x[3,b-1,i]-x[3,b,i])^2
        end
    end
    return A
end

function Path_X(A,x)
    B = size(x,2)
    N = size(x,3)
    for j in 1:N
        for i in 1:j
            A[i,j] += (x[1,1,i]-x[1,B,j])^2
            A[i,j] += (x[2,1,i]-x[2,B,j])^2
            A[i,j] += (x[3,1,i]-x[3,B,j])^2
    end end
    return A
end




