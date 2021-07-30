using Distributed
addprocs(4)
@everywhere include("get_objective_difference_solenoid.jl")
get_JavgN(0.001)
