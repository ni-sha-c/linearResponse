using JLD
using SharedArrays
using Distributed
include("../examples/lorenz.jl")
function obj_fun(x,y,z)
		return 	(z-28.0)^2.0
end
function obj_fun_erg_avg(s)
	nSteps = 100000
	u = 2*pi*rand(3)
	J = 0.
	u_trj = step(u,s,nSteps)
	x, y, z = view(u_trj,1,:), view(u_trj,2,:), view(u_trj,3,:)
	J = sum(obj_fun.(x,y,z)/nSteps)
	return J
end
function get_Javg_vs_s()
	n_pts = 25
	n_rep = 16000
	s = [10.0, 28.0, 8/3]
	s_ind = LinRange(27.0,32.0,n_pts)
	J = zeros(n_pts)
	J_proc = SharedArray{Float64}(n_rep)
	for i = 1:n_pts
		s[2] = s_ind[i]
		J_proc .= 0.
		t = @distributed for n=1:n_rep
			J_proc[n] = obj_fun_erg_avg(s)/n_rep
		end
		wait(t)
		J[i] = sum(J_proc) 
	end
	save("../data/obj_erg_avg/lorenz/zm28sq.jld",
		 "s2", s_ind,
		"J", J)
end


