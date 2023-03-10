using DPP

x = randn(2,500) #some points in dim 2

#compute a kernel matrix for the points in x
L = [ exp(-norm(a-b)^2) for a in eachcol(x), b in eachcol(x) ]
dpp = EllEnsemble(L) #form a L-ensemble based on the L matrix
rescale!(dpp,50) #scale so that the expected size is 50
ind = sample(dpp) #a sample from the DPP (indices)

using Plots
scatter(x[1,:],x[2,:],marker_z = map((v) -> v ∈ ind, 1:size(x,2)),legend=:none,alpha=.75) #show the selected points in white