using JLD
using SharedArrays
using Distributed
include("../examples/baker.jl")
function obj_fun(x,y)
		return cos(4*y)	
end
function obj_fun_erg_avg(s)
	nSteps = 10000
	u = 2*pi*rand(2)
	J = 0.
	u_trj = step(u,s,nSteps-1)
	x, y = view(u_trj,1,:), view(u_trj,2,:)
	J = sum(obj_fun.(x,y)/nSteps)
	return J
end
function get_Javg_vs_s(ind)
	s = zeros(4)
	n_pts = 30
	n_rep = 16000
	s_ind = LinRange(0.,0.4,n_pts)
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
	save("../data/obj_erg_avg/cos4y_s$(ind)_sens.jld",
		 "s$ind", s_ind,
		"J", J)
end


