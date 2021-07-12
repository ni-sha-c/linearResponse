using Distributed
addprocs(16)
@everywhere include("../src/get_sens_solenoid.jl")
si = LinRange(1.,2.0,20)
s = zeros(3,size(si)[1])
s[1,:] .= si
s[3,:] .= 0.0
s[2,:] .= 4.0
get_sens(s)
