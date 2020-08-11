include("../examples/baker.jl")
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function unstable_sens(s,nSteps)
	u = 2*pi*rand(2)
	u_trj = step(u, s, nSteps)'
	x, y = view(u_trj,:,1), view(u_trj,:,2)
	dJds = 0.
	N = 10
	nSteps = nSteps + 1
	g = zeros(nSteps)
	J = cos.(4*y)
	dXudx1 = zeros(2,nSteps)
	Xu = [sin.(x) zeros(nSteps)]'
	for i = 1:nSteps-1
		x1, x2 = x[i], y[i]
		z = 2.0 + s[1]*cos(x1)
		dzdx = -s[1]*sin(x1)
		dXudx1[:,i] = [cos(x1)/z, 0]   
		g[i+1] = g[i]/z - dzdx/z*z 
	end
	Xu = Xu[1,:]
	dXudx1 = dXudx1[1,:]
	for n = 1:N
		J_shift = J[n+1:end]
		nJ = length(J_shift)
		g_shift = g[2:nJ+1]
		dXudx1_shift = dXudx1[1:nJ]
		Xu_shift = Xu[1:nJ]
		dJds -= dot(J_shift,dXudx1_shift)/nJ 
		dJds -= dot(J_shift,g_shift.*Xu_shift)/nJ
	end
	#@show le, dJds
	return dJds
	
end
function get_unstable_sens(s1)
	p = 4
	nSteps = 5000
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
	save("../data/unstable_sens/dJds1.jld", "s1",
	     s1, "dJds", dJds)
end
	
