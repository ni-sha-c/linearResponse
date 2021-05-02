using Distributed
addprocs(16)
@everywhere include("get_objective_lorenz.jl")
get_Javg_vs_s()
