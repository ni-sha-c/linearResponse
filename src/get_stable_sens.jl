include("../examples/baker.jl")
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function stable_sens(s,nSteps)
	vs = zeros(2)
	q = rand(2)
	q /= norm(q)
	le = 0.
	u = 2*pi*rand(2)
	u_trj = step(u, s, nSteps-1)
	du_trj = dstep(u_trj, s)
	x = pert(u_trj, 4)
	dJds = 0.
	for i = 1:nSteps
		vs .= du_trj[:,:,i]*vs + x[:,i]
		q .= du_trj[:,:,i]*q 
		nm_q = norm(q)
		le += log(nm_q)/nSteps
		q ./= nm_q
		vs .-= dot(vs,q)*q
		x2 = (i==nSteps) ? step(u_trj[:,i],
			s, 1)[2,end] : u_trj[2,i+1]
		dJdu = [0., -4*sin(4*x2)]
		dJds += dot(dJdu, vs)/nSteps
	end
	#@show le, dJds
	return dJds
end
function get_stable_sens(s4)
	p = 4
	nSteps = 100
	s = zeros(p)
	# J = cos(4y)
	n_exps = size(s4)[1]
	dJds = zeros(n_exps)
	n_rep = 160000
	dJds_proc = SharedArray{Float64}(n_rep)
	for k=1:n_exps
		s[4] = s4[k]
		dJds_proc .= 0.
		t = @distributed for j=1:n_rep
			dJds_proc[j] = stable_sens(s,nSteps)/n_rep
		end
		wait(t)
		dJds[k] = sum(dJds_proc)
		@show dJds[k]
	end
	save("../data/stable_sens/dJds4.jld", "s4",
	     s4, "dJds", dJds)
end
	
