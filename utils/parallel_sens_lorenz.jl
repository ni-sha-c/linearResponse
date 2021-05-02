using Distributed
addprocs(16)
@everywhere include("../src/get_sens_lorenz.jl")
si = LinRange(27.5,30.5,10)
s = zeros(3,size(si)[1])
s[2,:] .= si
s[1,:] .= 10.
s[3,:] .= 8.0/3
get_sens(s)
