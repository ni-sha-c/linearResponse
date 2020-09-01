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
	@show size(pp)
	dJds_st = 0.
	dJds_ust = 0.
	N = 4
	nSteps = nSteps + 1
	g = zeros(nSteps)
	J = cos.(4*y)
	
	# Df represents the derivative of f 
	# on the 
	# unstable manifold.
	# Xu = a q
	Xu = zeros(2, nSteps)
	a = zeros(nSteps)
	Da = zeros(nSteps)
	Dq = zeros(2)
	Dvs = zeros(2)
	

	for i = 1:nSteps-1
		dui = du_trj[:,:,i]
		ppi = pp[:,i]
		xi, yi = x[i], y[i]
		dppi = [cos(xi) 0; cos(xi)*sin(yi) sin(xi)*cos(yi)]
		vs .= dui*vs + ppi
		q .= dui*q
		z = norm(q)
		z2 = z*z
		q ./= z
		a[i] = dot(vs, q)
		Xu[:, i] = a[i]*q 
		vs .-= dot(vs,q)*q
		d2q = reshape([dot(ddu_trj[:,j,i],
					q) for j=1:4],2,2)
		dz = d2q*q/z2 + dui*Dq/z2
		dzdx = dot(dz, q)
		Dq .= dz - dzdx*q
		Dvs .= d2q*vs/z + dui*vs/z + 
				dppi*q - Da[i]*q - a[i]*Dq
		Delta_Da = dot(Dvs, q) - dot(vs, Dq) 
		Dvs .+= Delta_Da*q
		Da[i] = dot(q, Dq) + dot(d2q*vs, q)/z + 
				dot(dui*vs/z, q) + dot(dppi*q, q) + 
				Delta_Da
		g[i+1] = g[i]/z - dzdx/z*z 
		dJdu = [0., -4*sin(4*yi)]
		dJds_st += dot(dJdu, vs)/nSteps
	end
	Xu = Xu[1,:]
	#dXuq = dXuq[1,:]
	for n = 1:N
		J_shift = J[n+1:end]
		nJ = length(J_shift)
		g_shift = g[2:nJ+1]
		Da_shift = Da[1:nJ]
		Xu_shift = Xu[1:nJ]
		dJds_ust -= dot(J_shift,Da_shift)/nJ 
		dJds_ust -= dot(J_shift,g_shift.*Xu_shift)/nJ
	end
	#@show le, dJds
	return dJds_st + dJds_ust
end

function get_sens(s)
	nSteps = 5000
	# J = cos(4y)
	n_exps = size(s1)[1]
	dJds = zeros(n_exps)
	n_rep = 1600
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
	     s1, "dJds", dJds)
end
	
