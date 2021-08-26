include("../examples/solenoid.jl")
using LinearAlgebra
using JLD
using Zygote
function sens(s,nSteps)
    d = 3
    s_ind = 3
    u = rand(d)
    
    n_runup = 500
    for i = 1:n_runup
        u = next(u,s)
    end
    
    dui = rand(d,d)
    ppi = rand(d)
    dppi = rand(d,d)

    dpert(x) = jacobian(x -> pert(x,s,s_ind), x)[1]
    d2next(x) = jacobian(x -> dstep(x,s), x)[1]
    
    vs = zeros(d)
    q = rand(d)
    q /= norm(q)

    pp = rand(d)
    q1 = rand(d)
    vs1 = zeros(d)
    
    dJds_st = 0.
    dJds_ust = 0.
    N = 12
    nSteps = nSteps + 1
    g = zeros(nSteps)
    
    
    # Df represents the derivative of f 
    # on the 
    # unstable manifold.
    # Xu = a q
    # d = 2
    Xu = zeros(d, nSteps)
    a = zeros(nSteps)
    J = zeros(nSteps)
    Da = zeros(nSteps)
    Dq = zeros(d)
    Dvs = zeros(d)
    Dvs1 = zeros(d)
    normvs = zeros(nSteps)

    # nSteps large number.
    for i = 1:nSteps-1
        J[i] = u[1]*u[1] + u[2]*u[2]    
        dJdu = [2*u[1], 2*u[2], 0.]
        ppi .= pert(u,s,s_ind)
        dui .= dstep(u,s) # for large systems, can't store Jacobian.
        dppi .= dpert(u)*inv(dui) # profile against using solve (backslash) here.                   
        q1 .= dui*q
        z = norm(q1)
        q1 ./= z

        z2 = z*z
        
        vs1 .= dui*vs + ppi
        a[i+1] = dot(vs1, q1)
        vs1 .-= a[i+1]*q1

        ddui = d2next(u)        
        d2q = reshape(ddui*q,d,d)
        Dq = d2q*q/z2 + dui*Dq/z2
        dzdx = dot(z*Dq, q1)
        Dq .= Dq .- dot(Dq,q1)*q1
        


        Dvs1 .= d2q*vs/z + dui*Dvs/z + dppi*q1
        Da[i+1] = dot(vs1, Dq) + dot(Dvs1, q1)
        Dvs1 .= Dvs1 - Da[i+1]*q1 - a[i+1]*Dq
        


        g[i+1] = g[i]/z - dzdx/z 

        dJds_st += dot(dJdu, vs)/nSteps
        
		normvs[i+1] = norm(vs1)

        vs .= vs1
        q .= q1
        Dvs .= Dvs1
        u .= next(u,s)

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
    save("../data/sens/solenoid/metrics_$(nSteps).jld", 
		 "a", a, "normvs", normvs, "b", Da, "g", g, 
		 "dJds_st", dJds_st, "dJds_ust", dJds_ust)
	return dJds_st + dJds_ust
end

