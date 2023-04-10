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

struct Potential{Float64}
    _PES
    Dimension::Int
    scale_begin::Vector{Float64}
    scale_end::Vector{Float64}
    bin::Vector{Float64}
end


