using Distributed
addprocs(16)
@everywhere include("../src/get_sens_solenoid.jl")
si = LinRange(0.,0.,1)
s = zeros(3,size(si)[1])
s[1,:] .= 1.0
s[3,:] .= si
s[2,:] .= 4.0
get_sens(s)
