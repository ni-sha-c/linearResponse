dt = 0.005
function lorenz63(u0, s, n)
    sigma, rho, beta = s
    d, m = size(u0)
    n = n+1
    u_trj = zeros((m,d,n))
    u_trj[:,:,1] = u0'
    for i = 2:n
    	x = u_trj[:,1,i-1]
    	y = u_trj[:,2,i-1]
    	z = u_trj[:,3,i-1]
    
    	u_trj[:,1,i] = x + dt*(sigma.*(y - x))
    	u_trj[:,2,i] = y + dt*(x.*(rho .- z) - y)
    	u_trj[:,3,i] = z + dt*(x.*y - beta.*z)
    end 
    return permutedims(u_trj,[3,2,1])
end
function dlorenz63(u, s)
    sigma, rho, beta = s
    n, d = size(u)
    x = view(u,:,1)
    y = view(u,:,2)
    z = view(u,:,3)
    du = zeros(n, d, d)
    @. du[:,1,1] = 1.0 - dt*sigma
    @. du[:,1,2] = dt*sigma
    @. du[:,2,1] = dt*(rho - z) 
    @. du[:,2,2] = 1.0 - dt
    @. du[:,2,3] = -dt*x 
    @. du[:,3,1] = dt*y
    @. du[:,3,2] = dt*x
    @. du[:,3,3] = 1.0 - dt*beta
    return reshape([du[:,:,1]'; du[:,:,2]'; 
    				du[:,:,3]'], d, d, n)
end
function perturbation(u,s)
    n, d = size(u)
    # the perturbation in row i in T_{u_(i+1)} M
    return [zeros(1,n); dt*u[:,1]'; zeros(1,n)]
end
function vectorField(u,s)
    n, d = size(u)
    sigma, rho, beta = s
    u = u'
    x, y, z = u[1,:], u[2,:], u[3,:]
    return [sigma.*(y - x)  x.*(rho .- z) - y  x.*y - beta.*z]'
end
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
function step(x0, s, n)
	x_trj = zeros(3, n+1)
	x_trj[:,1] = x0
	sigma, rho, beta = s
    d, m = size(u0)
    n = n+1
    for i = 2:n
    	x = x_trj[1,i-1]
    	y = x_trj[2,i-1]
    	z = x_trj[3,i-1]
    
    	x_trj[1,i] = x + dt*(sigma*(y - x))
    	x_trj[2,i] = y + dt*(x*(rho - z) - y)
    	x_trj[3,i] = z + dt*(x*y - beta*z)
    end 
    return x_trj
end
function next(u, s)
	x, y, z = u[1], u[2], u[3]
	u_next = similar(u)
	sigma, rho, beta = s
	u_next[1] = x + dt*(sigma*(y - x))
    u_next[2] = y + dt*(x*(rho - z) - y)
    u_next[3] = z + dt*(x*y - beta*z)
	return u_next
end

function dstep(u::Array{Float64,1},s::Array{Float64,1})
	du = zeros(3,3)
	x, y, z = u[1], u[2], u[3]
	sigma, rho, beta = s
    @. du[1,1] = 1.0 - dt*sigma
    @. du[1,2] = dt*sigma
    @. du[2,1] = dt*(rho - z) 
    @. du[2,2] = 1.0 - dt
    @. du[2,3] = -dt*x 
    @. du[3,1] = dt*y
    @. du[3,2] = dt*x
    @. du[3,3] = 1.0 - dt*beta
	return du
end

function dstep(u::Array{Float64,2},s::Array{Float64,1})
	n = size(u)[2]
	du = zeros(3,3,n)
	x, y, z = view(u,1,:),view(u,2,:),view(u,3,:)
	sigma, rho, beta = s
	@. du[1,1,:] = 1.0 - dt*sigma
    @. du[1,2,:] = dt*sigma
    @. du[2,1,:] = dt*(rho - z) 
    @. du[2,2,:] = 1.0 - dt
    @. du[2,3,:] = -dt*x 
    @. du[3,1,:] = dt*y
    @. du[3,2,:] = dt*x
    @. du[3,3,:] = 1.0 - dt*beta
	return du
end
function pert(u::Array{Float64,2}, p::Int64)
	n = size(u)[2]
	x, y = view(u,1,:), view(u,2,:)
	sx, sy = sin.(x), sin.(2*y)./2
	if p==1
		return reshape([sx zeros(n)],n,2)'
	elseif p==2
		return reshape([sy.*sx zeros(n)],n,2)' 
	elseif p==3
		return reshape([zeros(n) sy.*sx],n,2)' 
	elseif p==4
		return reshape([zeros(n) sy],n,2)'
	else
		println("parameter perturbation is only defined 
				for indices 1 through 4")
		return 0
	end
end
function pert(u::Array{Float64,1}, p::Int64)
	x, y = u[1], u[2]
	sx, sy = sin(x), sin(2*y)/2
	if p==1
		return [sx, 0]
	elseif p==2
		return [sy.*sx, 0]
	elseif p==3
		return [0, sy.*sx] 
	elseif p==4
		return [0, sy]
	else
		println("parameter perturbation is only defined 
				for indices 1 through 4")
		return 0
	end
end
function dpert_x(u::Array{Float64,1},s::Array{Float64,1})
	# here we asssume s = [a, a, 0, a]
	x, y = u[1], u[2]
	sx, sy = sin(x), sin(2*y)/2
	cx = cos(x)
	dxu = 2 .+ s[1]*cx .+ s[2]*sy.*cx 
	dxpx1 = cx + cx*sy 
	return [dxpx1/dxu, 0.]
end
function d2step(u::Array{Float64,2}, s::Array{Float64,1})
	n = size(u)[2]
	ddu = zeros(2,4,n)
	for i = 1:n
		x, y = u[1,i], u[2,i]
		sx, sy = sin(x), sin(2*y)/2
		dsx, dsy = cos(x), cos(2*y)
		dxdsx, dydsy = -sin(x), -2*sin(2*y)
		ddu[:,1,i] = [s[1]*dxdsx + s[2]*sy*dxdsx, s[2]*dsy*dsx] 
		ddu[:,2,i] = [s[3]*dxdsx*sy, s[3]*dsx*dsy]
		ddu[:,3,i] = [s[2]*dsx*dsy, s[2]*sx*dydsy]
		ddu[:,4,i] = [s[3]*dsx*dsy, s[4]*dydsy + s[3]*sx*dydsy]
	end
	return ddu
end
function d2step(u::Array{Float64,1}, s::Array{Float64,1})
	ddu = zeros(2,4)
	x, y = u[1], u[2]
	sx, sy = sin(x), sin(2*y)/2
	dsx, dsy = cos(x), cos(2*y)
	dxdsx, dydsy = -sin(x), -2*sin(2*y)
	ddu[:,1] = [s[1]*dxdsx + s[2]*sy*dxdsx, s[2]*dsy*dsx] 
	ddu[:,2] = [s[3]*dxdsx*sy, s[3]*dsx*dsy]
	ddu[:,3] = [s[2]*dsx*dsy, s[2]*sx*dydsy]
	ddu[:,4] = [s[3]*dsx*dsy, s[4]*dydsy + s[3]*sx*dydsy]
	return ddu
end
