using AdvancedMH
N = 2
B = 10

lp(x) = begin 
    p = reshape(x,N,B)
    return U(p) + fermi(p) 
end

function AMC(P,steps::Real,lp)

    density_model = AdvancedMH.DensityModel(x -> lp(x))
    # proposal = Array{Normal{Float64}}(undef, N,dim,beads)

    # for i in 1:N
    #     proposal[i,:,:] .= Normal(0, 1)
    # end
    proposal = RandomWalkProposal(MvNormal(zeros(N,dim,beads)[:], I))
    sampler = AdvancedMH.MetropolisHastings(proposal)

    chain = AdvancedMH.sample(density_model, sampler, convert(Int, steps);
                            init_params=P[:])

    return chain
end