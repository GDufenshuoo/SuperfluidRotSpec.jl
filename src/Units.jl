
"""
This code defines a struct Unit which has several fields of type F, where F is a subtype of Real. 
    The fields represent various physical constants and are named 
    according to their usual symbols.

The fields are:

`cᵁ`: the speed of light in a vacuum
`ħ`: the reduced Planck constant
`eᶜ`: the elementary charge
`mₑ`: the mass of an electron
`a₀`: the Bohr radius
`sᵁ`: the Stefan-Boltzmann constant
`Eᵁ`: the energy of the unit
`Eᵁₖ`: the Boltzmann constant
"""
struct Unit{F<:Real}
    cᵁ::F
    ħ::F
    eᶜ::F
    mₑ::F
    a₀::F
    sᵁ::F
    Eᵁ::F

    Eᵁₖ::F
end

"""
## In Unit of Atomic Unit.

The constant Atomicᵁ is an instance of the Unit struct with specific values for each field. 
"""
const Atomicᵁ = Unit(1.,1.,1.,1.,1.,1.,1.,3.1668105084779793e-6)


