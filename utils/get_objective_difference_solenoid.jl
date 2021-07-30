using JLD
using SharedArrays
using Distributed
include("../examples/solenoid.jl")
function obj_fun(u)
	return u[1]^2 + u[2]^2	
end
function obj_fun_erg_avg(s, nSteps)
	nRunup = 500
	u = rand(3)
	J = 0.
	for i = 1:nRunup
		u = next(u,s)
	end
	for i = 1:nSteps
		J += obj_fun(u)/nSteps
		u = next(u,s)
	end
	return J
end
function get_JavgN(sind)
	s = [1.0, 4.0, sind]
	n_rep = 160
	n_pts = 8
	J = zeros(n_pts)
	J_proc = SharedArray{Float64}(n_rep)
	J_proc .= 0.
	Ni = 1200
	N = zeros(Int64, n_pts)
	for i = 1:n_pts
		N[i] = Ni*n_rep
		@show Ni
		t = @distributed for n=1:n_rep
			J_proc[n] = obj_fun_erg_avg(s,Ni)/n_rep
		end
		wait(t)
		J[i] = sum(J_proc) 
		Ni = Ni*5
	end
	save("../data/finiteDifference/solenoid/r2_$(sind).jld","s", sind,"J", J, "N", N)
end
