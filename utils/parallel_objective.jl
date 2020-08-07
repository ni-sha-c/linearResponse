using Distributed
addprocs(16)
@everywhere include("get_objective.jl")
get_Javg_vs_s(4)
