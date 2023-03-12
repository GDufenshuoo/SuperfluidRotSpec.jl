using JLD2
using Interpolations

function set_potention(Pot)

    if Pot.Dimension == 1
        scale = (Pot.scale_begin[1]:Pot.bin[1]:Pot.scale_end[1])
    elseif Pot.Dimension == 2
        scale = (Pot.scale_begin[1]:Pot.bin[1]:Pot.scale_end[1], Pot.scale_begin[2]:Pot.bin[2]:Pot.scale_end[2])
    end

    return LinearInterpolation(
        scale,
        Pot._PES,
        extrapolation_bc=Line())
end

struct Potential_{Float64}
    _PES
    Dimension::Int
    scale_begin::Vector{Float64}
    scale_end::Vector{Float64}
    bin::Vector{Float64}
end

const am2An = 5.29177210903e-1
const pH2_H2 = set_potention(load("OCS.PES")["paraH2_paraH2"])
const OCS_pH2 = set_potention(load("OCS.PES")["OCS_paraH2"])

function 𝑈(x,N,B)
    U1 = 0.0
    U2 = 0.0
    @floop for i in 1:N
        for b in 1:B
            r = norm(x[:,b,i])*am2An
            cos = x[1,b,i]/r
            @reduce U1 += (r > 30.0 ? 0.0 : OCS_pH2(r,cos))
    end end
    @floop for i in 2:N
        for j in 1:i
            for b in 1:B
                @reduce U2 += pH2_H2(norm(x[:,b,i].-x[:,b,j])*am2An)
    end end end
    return (U1+U2)/B
end

