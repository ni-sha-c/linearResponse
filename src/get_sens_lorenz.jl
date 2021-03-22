include("../examples/lorenz.jl")
using LinearAlgebra
using JLD
using SharedArrays
using Distributed
function sens(s,nSteps)
    d = 3
    n_runup = 2000
    u = rand(d)
    for i = 1:n_runup
	u = next(u, s)
    end
    u_trj = step(u, s, nSteps)
    x, y, z = view(u_trj,1,:), view(u_trj,2,:), view(u_trj,3,:)
   
    du_trj = dstep(u_trj, s)
    ddu_trj = d2step(u_trj,s)
    vs = zeros(d)
    q = rand(d)
    q /= norm(q)

    pp = pert(u_trj, s)
    dpp = dpert(u_trj, s)
    q1 = rand(d)
    vs1 = zeros(d)
    
    dJds_st = 0.
    dJds_ust = 0.
    N = 12
    nSteps = nSteps + 1
    g = zeros(nSteps)
    J = copy(z)
   
    # Df represents the derivative of f 
    # on the 
    # unstable manifold.
    # Xu = a q
    Xu = zeros(d, nSteps)
    a = zeros(nSteps)
    Da = zeros(nSteps)
    Dq = zeros(d)
    Dvs = zeros(d)
    Dvs1 = zeros(d)
    rho = s[2]
    # nSteps large number.
    for i = 1:nSteps-1
	@show i
	dui = du_trj[:,:,i] # for large systems, can't store Jacobian.
    	ppi = pp[:,i]
	dppi = dpp[:,:,i]
	xi, yi, zi = x[i], y[i], z[i]

    	q1 .= dui*q
    	alpha = norm(q1)
    	q1 ./= alpha
    	alpha2 = alpha*alpha
    	
    	vs1 .= dui*vs + ppi
    	a[i+1] = dot(vs1, q1)
    	vs1 .-= a[i+1]*q1
    	
    	d2q = reshape([dot(ddu_trj[:,j,i],
    				q) for j=1:d*d],d,d)
    	Dq = d2q*q/alpha2 + dui*Dq/alpha2
    	dalphadx = dot(alpha2*Dq, q1)
    	Dq .= Dq .- dot(Dq,q1)*q1

    	Dvs1 .= d2q*vs/alpha + dui*Dvs/alpha + dppi*q1/alpha
    	Da[i+1] = dot(vs1, Dq) + dot(Dvs1, q1)
    	Dvs1 .= Dvs1 - Da[i+1]*q1 - a[i+1]*Dq
    	Delta_Da = dot(Dvs1, q1) + dot(vs1, Dq) 
    	Dvs1 .-= Delta_Da*q
    	Da[i+1] += Delta_Da
    	
    	g[i+1] = g[i]/alpha - dalphadx/alpha2 
    	dJdu = [0., 0., 1.]
    	dJds_st += dot(dJdu, vs)/nSteps
    	vs .= vs1
    	q .= q1
    	Dvs .= Dvs1
    end

    for n = 1:N
    	J_shift = J[n+1:end]
    	nJ = length(J_shift)
    	g_shift = g[2:nJ+1]
    	Da_shift = Da[2:nJ+1]
    	a_shift = a[2:nJ+1]
    	dJds_ust -= dot(J_shift,Da_shift)/nJ 
    	dJds_ust -= dot(J_shift,g_shift.*a_shift)/nJ
    end
   
    @show sum(J)/nSteps
    return dJds_st + dJds_ust
end
function get_sens(s)
    nSteps = 500000
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
    save("../data/sens/lorenz/dJds.jld", "rho",
         s[2,:], "dJds", dJds)
end
    
