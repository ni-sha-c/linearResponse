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
    f_trj = flow(u_trj,s)

    vs = zeros(d)
    q = rand(d)
    q /= norm(q)

    pp = pert(u_trj, s)
    dpp = dpert(u_trj, s)
    q1 = rand(d)
    vs1 = zeros(d)
    
    dJds_st = 0.
    dJds_ust = 0.
    dJds_c = 0.
    N = 1000
    nSteps = nSteps + 1
    g = zeros(nSteps)
    J = (z .- 28.0).^2.0
    #dJdu = [0., 0., 1.]
    # Df represents the derivative of f 
    # on the 
    # unstable manifold.
    a = zeros(nSteps)
    c = zeros(nSteps)
    beta = zeros(nSteps)
    dJf = zeros(nSteps)
    Da = zeros(nSteps)
    Dc = zeros(nSteps)
    Dq = zeros(d)
    Dq1 = zeros(d)
    Dvs = zeros(d)
    Dvs1 = zeros(d)
    d2varphi = zeros(d,d,d)
    d2q = zeros(d,d)
    rho = s[2]
    
    # nSteps large number.
    for i = 1:nSteps-1
        dui = du_trj[:,:,i] # for large systems, can't store Jacobian.
        ppi = pp[:,i]
        dppi = dpp[:,:,i]
    	beta[i+1] = norm(f_trj[:,i+1])
    	f1 = f_trj[:,i+1]/beta[i+1]


        q1 .= dui*q
        alpha = norm(q1)
        q1 ./= alpha
        alpha2 = alpha*alpha
    	th = dot(q1, f1)
        th2 = 1.0 - th*th
        r1 = dunit_flow([x[i+1],y[i+1],z[i+1]],s)*q1

        vs1 .= dui*vs + ppi
        a[i+1] = dot(vs1, q1 - th*f1)/th2
    	c[i+1] = dot(vs1, f1 - th*q1)/th2

        d2varphi .= permutedims(ddu_trj[:,:,:,i], 
				 [1,3,2]) 
    	for j = 1:d
    		d2q[:,j] .= d2varphi[:,:,j]*q
    	end
         
        Dq1 .= d2q*q/alpha2 + dui*Dq/alpha2
        dalphadx = dot(Dq1, q1)*alpha
	Dq1 .= Dq1 .- dalphadx*q1/alpha
    	
        Dvs1 .= d2q*vs/alpha + dui*Dvs/alpha + dppi*q1
    	Dth = dot(q1, r1) + dot(f1, Dq1)
    	constant = 2*th*Dth/th2/th2
    	Da[i+1] = constant*dot(vs1, q1 - th*f1)
    	Da[i+1] += 1/th2*(dot(vs1, Dq1 - Dth*f1 - th*r1) + 
    				 dot(Dvs1, q1 - th*f1))

    	Dc1 = constant*dot(vs1, f1 - th*q1)
    	Dc1 += 1/th2*(dot(vs1, r1 - Dth*q1 - th*Dq1) + 
    				 dot(Dvs1, f1 - th*q1))

       
        g[i+1] = g[i]/alpha - dalphadx/alpha 
        vs1 .= vs1 .- c[i+1]*f1 .- a[i+1]*q1

    	Dvs1 .= Dvs1 .- a[i+1]*Dq1 .- Da[i+1]*q1 .- c[i+1]*r1 .- Dc1*f1
	
	dJdu = [0.,0., 2*(z[i] - 28.0)]
        dJds_st += dot(dJdu, vs)/nSteps
    	dJf[i+1] = dot(dJdu, f1)

        vs .= vs1
        q .= q1
        Dvs .= Dvs1
	Dq .= Dq1
    end

    for n = 1:N
        J_shift = J[n+1:end]
    	dJf_shift = dJf[n+1:end]
    	betan_shift = beta[n+1:end]
    	
        nJ = length(J_shift)
        g_shift = g[2:nJ+1]
        Da_shift = Da[2:nJ+1]
        a_shift = a[2:nJ+1]
    	c_shift = c[2:nJ+1]
    	beta_shift = beta[2:nJ+1]

        dJds_ust -= dot(J_shift,Da_shift)/nJ 
        dJds_ust -= dot(J_shift,g_shift.*a_shift)/nJ
    	dJds_c += dot(c_shift./beta_shift, dJf_shift.*betan_shift)/nJ 
    end
    return dJds_st + dJds_ust + dJds_c
    
end
function get_sens(s)
    nSteps = 1000000
    n_exps = size(s)[2]
    dJds = zeros(n_exps)
    var_dJds = zeros(n_exps)
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
        var_dJds[k] = n_rep*sum(dJds_proc.*dJds_proc) - (dJds[k])^2
    end
    save("../data/sens/lorenz/dJds_squareObj.jld", "rho",
         s[2,:], "dJds", dJds, "var_dJds", var_dJds)
end
    
