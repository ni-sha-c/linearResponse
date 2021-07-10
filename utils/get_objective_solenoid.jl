using JLD
using SharedArrays
using Distributed
include("../examples/solenoid.jl")
function obj_fun(u)
		return cos(4*u[1])	
end
function obj_fun_erg_avg(s)
	nSteps = 10000
	nRunup = 500
	u = rand(d)
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
function get_Javg_vs_s(ind)
	s = [1.0, 4.0, 0.1]
	n_pts = 10
	n_rep = 16000
	s_ind = LinRange(0.,1.0,n_pts)
	J = zeros(n_pts)
	J_proc = SharedArray{Float64}(n_rep)
	for i = 1:n_pts
		@show s_ind[i]
		s[ind] = s_ind[i]
		#s[3] = s_ind[i]
		J_proc .= 0.
		t = @distributed for n=1:n_rep
			J_proc[n] = obj_fun_erg_avg(s)/n_rep
		end
		wait(t)
		J[i] = sum(J_proc) 
	end
	save("../data/obj_erg_avg/solenoid/cos4x_s$(ind)_sens.jld",
		 "s$ind", s_ind,
		"J", J)
end


