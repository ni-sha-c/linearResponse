include("../examples/baker.jl")
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function unstable_sens(s,nSteps)
	q = rand(2)
	q /= norm(q)
	le = 0.
	u = 2*pi*rand(2)
	u_trj = step(u, s, nSteps-1)
	du_trj = dstep(u_trj, s)
	Xu = pert(u_trj, 4)
	dJds = 0.
	N = 10
	g = zeros(nSteps+1)
	J = zeros(nSteps+1)
	dXudx1 = zeros(nSteps+1)
	dXudx1[1:nSteps] = cos.(u_trj[1:,])
	dXudx1[end] = cos(step(u_trj[:,end],s,1)[1,end])
	for i = 1:nSteps
		x1, x2 = (i==nSteps) ? step(u_trj[:,i],
			s, 1)[:,end] : u_trj[:,i+1]
		J[i+1] = cos(4*x2)
		z = 2.0 + s[1]*cos(x1)
		dzdx = -s[1]*sin(x1)
		g[i+1] = g[i]*z + dzdx 
	end

	for n = 1:N
		J_shift = J[n+1:end]
		nJ = length(J_shift)
		g_shift = g[2:nJ+1]
		dXudx1_shift = dXudx1[1:nJ]
		Xu_shift = Xu[1:nJ]
		dJds -= dot(J_shift,dXudx1_shift) 
		dJds -= dot(J_shift,g_shift.*Xu_shift)
	end
	#@show le, dJds
	return dJds
end
function get_unstable_sens(s1)
	p = 4
	nSteps = 100
	s = zeros(p)
	# J = cos(4y)
	n_exps = size(s1)[1]
	dJds = zeros(n_exps)
	n_rep = 160000
	dJds_proc = SharedArray{Float64}(n_rep)
	for k=1:n_exps
		s[1] = s1[k]
		dJds_proc .= 0.
		t = @distributed for j=1:n_rep
			dJds_proc[j] = unstable_sens(s,nSteps)/n_rep
		end
		wait(t)
		dJds[k] = sum(dJds_proc)
		@show dJds[k]
	end
	save("../data/stable_sens/dJds1.jld", "s1",
	     s1, "dJds", dJds1)
end
	
