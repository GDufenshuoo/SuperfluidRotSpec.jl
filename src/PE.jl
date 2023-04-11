
"""
am2An = 5.29177210903e-1
"""
function set_potention(Pot;U2u=1.0)

    if Pot.Dimension == 1
        scale = (Pot.scale_begin[1]*U2u:Pot.bin[1]*U2u:Pot.scale_end[1]*U2u)
    elseif Pot.Dimension == 2
        scale = (Pot.scale_begin[1]*U2u:Pot.bin[1]*U2u:Pot.scale_end[1]*U2u, 
        Pot.scale_begin[2]*U2u:Pot.bin[2]*U2u:Pot.scale_end[2]*U2u)
    end

    return LinearInterpolation(
        scale,
        Pot._PES,
        extrapolation_bc=Line())
end

struct Potential{Float64}
    _PES
    Dimension::Int
    scale_begin::Vector{Float64}
    scale_end::Vector{Float64}
    bin::Vector{Float64}
end
