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
	n_rep = 16
	n_pts = 2
	J = zeros(n_pts)
	J_proc = SharedArray{Float64}(n_rep)
	J_proc .= 0.
	Ni = 100
	for i = 1:n_pts
		t = @distributed for n=1:n_rep
			J_proc[n] = obj_fun_erg_avg(s,Ni)/(n_rep*Ni)
		end
		wait(t)
		J[i] = sum(J_proc) 
		Ni = Ni*5
	end
	save("../data/obj_erg_avg/solenoid/r2_$(sind).jld",
		 "s", s_ind,
		"J", J)
end


