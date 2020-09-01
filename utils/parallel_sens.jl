using Distributed
addprocs(16)
@everywhere include("../src/get_unstable_sens.jl")
si = LinRange(0.,1.,5)
s = zeros(4,5)
s[1,:] = si
s[3,:] = si
get_unstable_sens(s)
