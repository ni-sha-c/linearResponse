using Distributed
addprocs(16)
@everywhere include("plot_density.jl")
get_dist()
