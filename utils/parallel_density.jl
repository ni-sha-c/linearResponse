using Distributed
addprocs(16)
@everywhere include("get_srb.jl")
get_dist()
