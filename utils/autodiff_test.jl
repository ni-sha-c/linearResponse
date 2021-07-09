include("../examples/solenoid.jl")
using Zygote
using Test

s = [1.0,4.0,0.1]
x = rand(3)
g(x) = jacobian(x -> dstep(x,s), x)[1]
@show g(rand(3))

@time for i = 1:100
	x .= rand(3)
	du_ad = g(x)
end

@time for i = 1:100 
	x .= rand(3)
	du_ana = dstep(x, s)
end
#@test du_ad ≈ du_ana atol=1.e-8


