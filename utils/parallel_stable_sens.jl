using Distributed
addprocs(16)
@everywhere include("../src/get_stable_sens.jl")
s = LinRange(0.,1.,25)
get_stable_sens(s)
