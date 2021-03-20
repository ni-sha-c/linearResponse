using Distributed
addprocs(16)
@everywhere include("../src/get_sens.jl")
si = LinRange(0.,1.0,10)
s = zeros(4,10)
s[2,:] = si
#s[3,:] = si
get_sens(s)
