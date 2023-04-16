function Observe(lP,l::ClassicRotor)
    @unpack N = l
    B = 1
    P = size(lP,1)
    pb = zeros(2,B,N,P)
    for p in 1:P
        x = reshape(lP[p],3,B,N)
        for i in 1:N
            for b in 1:B
                r = norm(x[:,b,i])
                cos = x[1,b,i]/r
                pb[:,b,i,p] .= r*sqrt(1-cos^2),x[1,b,i]
    end end end
    return pb
end

function Observe(lP,l::SuperfluidFixRotor)
    @unpack N,B = l
    P = size(lP,1)
    pb = zeros(2,B,N,P)
    for p in 1:P
        x = reshape(lP[p],3,B,N)
        for i in 1:N
            for b in 1:B
                r = norm(x[:,b,i])
                cos = x[1,b,i]/r
                pb[:,b,i,p] .= r*sqrt(1-cos^2),x[1,b,i]
    end end end
    return pb
end

function Observe(lP,l::SuperfluidRotor)
    @unpack N,B = l
    P = size(lP,1)
    pb = zeros(2,B,N,P)
    φm = 3N*B
    for p in 1:P
        x = reshape(lP[p][begin:φm],3,B,N)
        for i in 1:N
            for b in 1:B
                r = norm(x[:,b,i])
                cos = x[1,b,i]/r
                pb[:,b,i,p] .= r*sqrt(1-cos^2),x[1,b,i]
    end end end
    return pb
end

function Observe(lP,l::SuperfluidRotor)
    @unpack N,B,rRB = l
    P = size(lP,1)
    pb = zeros(2,B,N,P)
    φm = 3N*B
    RB = fld(B,rRB)
    for p in 1:P
        x = reshape(lP[p][begin:φm],3,B,N)
        Rx = reshape(lP[p][φm+1:end],5,RB)
        for i in 1:N
            for b in 1:B
                rb_i = fld(b+rRB-1,rRB)
                rx = (x[:,b,i].-Rx[1:3,rb_i])
                r = norm(rx)
                cosθ = ix_rot_yz(rx,Rx[4:5,rb_i].%pi)/r

                rb_i = fld(b+rRB-1,rRB)
                rx = (x[:,b,i].-Rx[1:3,rb_i])
                r² = sum(abs2,rx)
                y = ix_rot_yz(rx,Rx[4:5,rb_i].%pi)
                pb[:,b,i,p] .= sqrt(r²-y^2),y
    end end end
    return pb
end


"""using Plots
# using Makie
histogram2d(pb[2,1,1,:],pb[1,1,1,:],bins=(40, 20), show_empty_bins=true,
normalize=:pdf,xlim = (-5, 5),ylim = (0, 5),background_color=:Black)

histogram2d(pb[2,1,2,:],pb[1,1,2,:],bins=(40, 20), show_empty_bins=true,
normalize=:pdf,xlim = (-5, 5),ylim = (0, 5),background_color=:Black)

pythonplot()"""
