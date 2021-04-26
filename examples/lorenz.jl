using LinearAlgebra
dt = 0.002
function flow(u::Array{Float64,2},s::Array{Float64,1})
    sigma, rho, beta = s
    x, y, z = view(u,1,:), view(u,2,:), view(u,3,:)
    return [sigma.*(y - x)  x.*(rho .- z) - y  x.*y - beta.*z]'
end
function flow(u::Array{Float64,1},s::Array{Float64,1})
    sigma, rho, beta = s
    x, y, z = u[1], u[2], u[3]
    return [sigma*(y - x), x*(rho - z) - y, x*y - beta*z]
end
function mag_flow(u::Array{Float64,2},s::Array{Float64,1})
	d, n = size(u)
	mags = zeros(n)
	for i = 1:n
		mags[i] = norm(flow(u[:,i],s))
	end
	return mags
end

function mag_flow(u::Array{Float64,1},s::Array{Float64,1})
	return norm(flow(u,s))
end
function dmag_flow(u::Array{Float64,1},s::Array{Float64,1})
	sigma, rho, beta = s
	x, y, z = u[1], u[2], u[3]
	b = mag_flow(u,s)
	dxb2 = -sigma*(y-x) + (x*(rho - z) - y)*(rho - z) + 
	        (x*y - beta*z)*y
	dyb2 = sigma*(y-x) - (x*(rho - z) - y) + 
	        (x*y - beta*z)*x
	dzb2 = -(x*(rho - z) - y)*x - 
			(x*y - beta*z)*beta
	return 1/b*[dxb2, dyb2, dzb2]
end
function step(u0, s, n)
    u_trj = zeros(3, n+1)
    u_trj[:,1] = u0
    n = n+1
    for i = 2:n
        u = u_trj[:,i-1]
    	k1 = dt*flow(u, s)
    	k2 = dt*flow(u .+ k1, s)
    	u_trj[:,i] = @. u + (k1 + k2)/2 	
    end 
    return u_trj
end
function next(u, s)
    k1 = dt*flow(u, s)
    k2 = dt*flow(u .+ k1, s)
    u_next = @. u + (k1 + k2)/2 
    return u_next
end
function dflow(u::Array{Float64,2}, s::Array{Float64,1})
    sigma, rho, beta = s
    d, n = size(u)
    x = view(u,1,:)
    y = view(u,2,:)
    z = view(u,3,:)
    du = zeros(n, d, d)
    @. du[:,1,1] = -sigma
    @. du[:,1,2] = sigma
    @. du[:,2,1] = (rho - z) 
    @. du[:,2,2] = -1.0 
    @. du[:,2,3] = -x 
    @. du[:,3,1] = y
    @. du[:,3,2] = x
    @. du[:,3,3] = -beta
    return reshape([du[:,:,1]'; du[:,:,2]'; 
            		du[:,:,3]'], d, d, n)
end
function dflow(u::Array{Float64,1},s::Array{Float64,1})
    du = zeros(3,3)
    x, y, z = u[1], u[2], u[3]
    sigma, rho, beta = s
    du[1,1] = -sigma
    du[1,2] = sigma
    du[2,1] = (rho - z) 
    du[2,2] = -1.0
    du[2,3] = -x 
    du[3,1] = y
    du[3,2] = x
    du[3,3] = -beta
    return du
end
function dstep(u::Array{Float64,2}, s::Array{Float64,1})
    d, n = size(u)
    du = zeros(d, d, n)
    dui = zeros(d,d)
    for i = 1:n
    	ui = u[:,i]
    	k1 = dt*flow(ui, s)
    	dk1 = dt*dflow(ui, s)
    	dk2 = dt*dflow(ui .+ k1, s)
    	dui .= 1.0*I(d) .+ 0.5*(dk1 .+ dk2*(1.0*I(d) .+ dk1))
    	du[:,:,i] .= dui
    end
    return du
end
function dstep(u::Array{Float64,1}, s::Array{Float64,1})
    d = 3
    du = zeros(d, d)
    k1 = dt*flow(u, s)
    dk1 = dt*dflow(u, s)
    dk2 = dt*dflow(u .+ k1, s)
    du = 1.0*I(d) .+ 0.5*(dk1 .+ dk2*(1.0*I(d) .+ dk1))
    return du
end

function d2flow(u::Array{Float64,2}, s::Array{Float64,1})
    n = size(u)[2]
    ddu = zeros(3,3,3,n)
    for i = 1:n
        ddu[2,3,1,i] = -1.0
    	ddu[3,2,1,i] = 1.0
        ddu[3,1,2,i] = 1.0
        ddu[2,1,3,i] = -1.0
    end
    return ddu
end
function d2flow(u::Array{Float64,1}, s::Array{Float64,1})
    ddu = zeros(3,3,3)
    ddu[2,3,1] = -1.0
    ddu[3,2,1] = 1.0
    ddu[3,1,2] = 1.0
    ddu[2,1,3] = -1.0
    return ddu
end
function d2step(u::Array{Float64,2}, s::Array{Float64,1})
    d, n = size(u)
    ddu = zeros(d, d, d, n)
    k1 = dt*flow(u, s)
    dk1 = dt*dflow(u,s)
    dk2 = dt*dflow(u .+ k1, s)
    d2 = dt*d2flow(u, s)
    d2p = permutedims(d2,[1, 3, 2, 4]) 
    A = zeros(d,d,d)
    for k = 1:n
    	for j = 1:d
    		A[:,:,j] .= d2p[:,:,j,k]*(1.0*I(d) .+ dk1[:,:,k]) 
    	end
    	A = permutedims(A, [1,3,2])
    	for i = 1:d
    			ddu[:,:,i,k] = 0.5*(d2[:,:,i,k] .+ 
    					dk2[:,:,k]*d2[:,:,i,k] .+ A[:,:,i]*
    					(1.0*I(d) .+ dk1[:,:,k]))	
    	end
    end
    return ddu
end
function d2step(u::Array{Float64,1}, s::Array{Float64,1})
    ddu = zeros(3, 3, 3)
    k1 = dt*flow(u, s)
    dk1 = dt*dflow(u,s)
    dk2 = dt*dflow(u .+ k1, s)
    d2 = dt*d2flow(u, s)
    d2p = permutedims(d2,[1, 3, 2]) 
    A = zeros(3,3,3)
    for j = 1:3
    	A[:,:,j] .= d2p[:,:,j]*(1.0*I(3) .+ dk1) 
    end
    A = permutedims(A, [1,3,2])
    for i = 1:3
    	ddu[:,:,i] = 0.5*(d2[:,:,i] .+ dk2*d2[:,:,i] .+ A[:,:,i]*
    						  (1.0*I(3) .+ dk1))	
    end
    
    return ddu
end
function pert(u::Array{Float64,2},s::Array{Float64,1})
    n = size(u)[2]
    dsflow(x) = [0., x, 0.]
	pp = zeros(3, n)
    for i = 1:n
		ui = u[:,i]
		k1 = dt*flow(ui, s)
		dk2 = dt*dflow(ui .+ k1,s)
		dsk1 = dt*dsflow(ui[1])
		dsk2 = dt*dsflow(ui[1] .+ k1[1]) .+ 
			   dk2*dsk1
		pp[:,i] = 0.5*(dsk1 + dsk2)
	end
	return pp
end
function pert(u::Array{Float64,1},s::Array{Float64,1})
	dsflow(x) = [0., x, 0.]
	k1 = dt*flow(u, s)

	dk2 = dt*dflow(u .+ k1,s)
	dsk1 = dt*dsflow(u[1])
	dsk2 = dt*dsflow(u[1] .+ k1[1]) .+ 
			   dk2*dsk1
	pp = 0.5*(dsk1 + dsk2)
	return pp
end
function dpert_next(u_arr::Array{Float64,2},s::Array{Float64,1})
	d, n = size(u_arr) 
	dpp = zeros(3,3,n)
	ddsflow = zeros(3,3)
	ddsflow[2,1] = 1.0
	dsflow(x) = [0., x, 0.]
	for j = 1:n
		u = u_arr[:,j]
		k1 = dt*flow(u, s)
		dk1 = dt*dflow(u,s)
		dk2 = dt*dflow(u .+ k1,s)
		ddk2 = dt*d2flow(u .+ k1, s)
		dsk1 = dt*dsflow(u[1])
		ddsk1 = dt*ddsflow
		d2 = dt*d2flow(u, s)
    	d2p = permutedims(d2,[1, 3, 2]) 
    	A = zeros(3,3,3)
    	for j = 1:3
    		A[:,:,j] .= d2p[:,:,j]*(1.0*I(3) .+ dk1) 
    	end
    	A = permutedims(A, [1,3,2])
		ddsk2 = dt*ddsflow*(1.0*I(3) .+ dk1) .+
			dk2*ddsk1 
		for i = 1:3
			ddsk2[:,i] .+= A[:,:,i]*dsk1
		end
		dpp[:,:,j] .= 0.5*(ddsk1 + ddsk2)			

  	end
	return dpp
end
function dpert_next(u::Array{Float64,1},s::Array{Float64,1})
    dpp = zeros(3,3)
	ddsflow = zeros(3,3)
	ddsflow[2,1] = 1.0
	
	k1 = dt*flow(u, s)
	dk1 = dt*dflow(u,s)
	
	dk2 = dt*dflow(u .+ k1,s)
	ddk2 = dt*d2flow(u .+ k1, s)
	
	dsflow(x) = [0., x, 0.]
	dsk1 = dt*dsflow(u[1])
	ddsk1 = dt*ddsflow
		
	d2 = dt*d2flow(u, s)
    d2p = permutedims(d2,[1, 3, 2]) 
    A = zeros(3,3,3)
    for j = 1:3
    	A[:,:,j] .= d2p[:,:,j]*(1.0*I(3) .+ dk1) 
    end
    A = permutedims(A, [1,3,2])

	ddsk2 = dt*ddsflow*(1.0*I(3) .+ dk1) .+
			dk2*ddsk1 
	for j = 1:3
		ddsk2[:,j] .+= A[:,:,j]*dsk1
	end
	dpp .= 0.5*(ddsk1 + ddsk2)			
	return dpp
end
function dpert(u::Array{Float64, 1}, s::Array{Float64,1})
	du = dstep(u,s)
	dpp_next = dpert_next(u,s)
	dpp = du'\dpp_next'
	return dpp'
end
function dpert(u::Array{Float64, 2}, s::Array{Float64,1})
	du = dstep(u,s)
	dpp_next = dpert_next(u,s)
	d, n = size(u)
	dpp = zeros(3,3,n)
	for i = 1:n
		dpp[:,:,i] .= du[:,:,i]'\dpp_next[:,:,i]'
		dpp[:,:,i] .= dpp[:,:,i]'
	end
	return dpp
end

### AD functions
function lorenz63_rhs_ad(du, u, s, t)
    du[1] = s[1]*(u[2] - u[1])
    du[2] = u[1]*(s[2] - u[3]) - u[2]
    du[3] = u[1]*u[2] - s[3]*u[3]
end
function lorenz63_ad(u0, s, n)
    t = n*dt
    prob = ODEProblem(lorenz63_rhs_ad, u0, (0.,t), s)
    sol = Array(solve(prob, Tsit5(), saveat=dt))
    return sol[:,end]
end
function obj_fun(u0, s)
    prob = ODEProblem(lorenz63_rhs_ad, u0, (0.,1.1), s)
    #_prob = remake(prob,u0=u0,p=s) 
    sol = solve(prob, Tsit5(), saveat=0.005)
    sum(sol[3,:])/size(sol,2)
end



