using Distributed
addprocs(16)
@everywhere include("../src/get_unstable_sens.jl")
s = LinRange(0.,1.,25)
get_unstable_sens(s)
