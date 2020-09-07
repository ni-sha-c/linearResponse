include("../examples/baker.jl")
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function sens(s,nSteps)
	u = 2*pi*rand(2)
	u_trj = step(u, s, nSteps)
	x, y = view(u_trj,1,:), view(u_trj,2,:)
	du_trj = dstep(u_trj, s)
	ddu_trj = d2step(u_trj,s)
	vs = zeros(2)
	q = rand(2)
	q /= norm(q)

	pp = pert(u_trj, 1) .+ pert(u_trj,3)
	q1 = rand(2)
	vs1 = zeros(2)
	
	dJds_st = 0.
	dJds_ust = 0.
	N = 11
	nSteps = nSteps + 1
	g = zeros(nSteps)
	J = cos.(4*y)
	
	# Df represents the derivative of f 
	# on the 
	# unstable manifold.
	# Xu = a q
	# d = 2
	Xu = zeros(2, nSteps)
	a = zeros(nSteps)
	Da = zeros(nSteps)
	Dq = zeros(2)
	Dvs = zeros(2)
	Dvs1 = zeros(2)
	
	# nSteps large number.
	for i = 1:nSteps-1
		dui = du_trj[:,:,i] # for large systems, can't store Jacobian.
		ppi = pp[:,i]
		xi, yi = x[i], y[i]
		dppi = [cos(xi) 0; cos(xi)*sin(yi) sin(xi)*cos(yi)]
		q1 .= dui*q
		z = norm(q1)
		q1 ./= z

		z2 = z*z
		
		vs1 .= dui*vs + ppi
		a[i+1] = dot(vs1, q1)
		vs1 .-= a[i+1]*q1

		
		d2q = reshape([dot(ddu_trj[:,j,i],
					q) for j=1:4],2,2)
		dz = d2q*q/z2 + dui*Dq/z2
		dzdx = dot(dz, q)
		Dq .= dz - dzdx*q
		
		Dvs1 .= d2q*vs/z + dui*Dvs/z + dppi*q1
		Da[i+1] = dot(vs1, Dq) + dot(Dvs1, q1)
		Dvs1 .= Dvs1 - Da[i+1]*q1 - a[i+1]*Dq
		Delta_Da = dot(Dvs1, q1) + dot(vs1, Dq) 
		Dvs .-= Delta_Da*q
		Da[i+1] += Delta_Da
		g[i+1] = g[i]/z - dzdx/z*z 
		dJdu = [0., -4*sin(4*yi)]
		dJds_st += dot(dJdu, vs)/nSteps

		vs .= vs1
		q .= q1
		Dvs .= Dvs1


	end
	for n = 1:N
		J_shift = J[n+1:end]
		nJ = length(J_shift)
		g_shift = g[2:nJ+1]
		Da_shift = Da[1:nJ]
		a_shift = a[1:nJ]
		dJds_ust -= dot(J_shift,Da_shift)/nJ 
		dJds_ust -= dot(J_shift,g_shift.*a_shift)/nJ
	end
	#@show le, dJds
	return dJds_st + dJds_ust
end
function get_sens(s)
	nSteps = 100000
	# J = cos(4y)
	n_exps = size(s)[2]
	dJds = zeros(n_exps)
	n_rep = 16
	dJds_proc = SharedArray{Float64}(n_rep)
	for k=1:n_exps
		sk = s[:,k]
		dJds_proc .= 0.
		t = @distributed for j=1:n_rep
			dJds_proc[j] = sens(sk,nSteps)/n_rep
		end
		wait(t)
		dJds[k] = sum(dJds_proc)
		@show dJds[k]
	end
	save("../data/sens/dJds.jld", "s1",
	     s[1,:], "dJds", dJds)
end
	
