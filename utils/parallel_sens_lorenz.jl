using Distributed
addprocs(16)
@everywhere include("../src/get_sens_lorenz.jl")
si = LinRange(25.,29.0,10)
s = zeros(3,10)
s[2,:] .= si
s[1,:] .= 10.
s[3,:] .= 8.0/3
get_sens(s)
