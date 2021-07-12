using Distributed
addprocs(16)
@everywhere include("get_objective_solenoid.jl")
get_Javg_vs_s(1)
