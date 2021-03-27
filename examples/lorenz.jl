dt = 0.001
function vectorField(u::Array{Float64,2},s::Array{Float64,1})
    sigma, rho, beta = s
    x, y, z = view(u,1,:), view(u,2,:), view(u,3,:)
    return [sigma.*(y - x)  x.*(rho .- z) - y  x.*y - beta.*z]'
end
function vectorField(u::Array{Float64,1},s::Array{Float64,1})
    sigma, rho, beta = s
	x, y, z = u[1], u[2], u[3]
    return [sigma.*(y - x)  x.*(rho .- z) - y  x.*y - beta.*z]'
end
function step(x0, s, n)
    x_trj = zeros(3, n+1)
    x_trj[:,1] = x0
    n = n+1
    for i = 2:n
        x = x_trj[:,i-1]
		k1 = dt*vectorField(x, s)
		k2 = dt*vectorField(x .+ k1, s)
    	@. x_trj[:,i] = x + (k1 + k2)/2 	
	end 
    return x_trj
end
function next(u, s)
    x, y, z = u[1], u[2], u[3]
    u_next = similar(u)
    return u_next
end
function dstep(u::Array{Float64,2}, s::Array{Float64,1})
    sigma, rho, beta = s
    d, n = size(u)
    x = view(u,1,:)
    y = view(u,2,:)
    z = view(u,3,:)
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
function dstep(u::Array{Float64,1},s::Array{Float64,1})
    du = zeros(3,3)
    x, y, z = u[1], u[2], u[3]
    sigma, rho, beta = s
    du[1,1] = 1.0 - dt*sigma
    du[1,2] = dt*sigma
    du[2,1] = dt*(rho - z) 
    du[2,2] = 1.0 - dt
    du[2,3] = -dt*x 
    du[3,1] = dt*y
    du[3,2] = dt*x
    du[3,3] = 1.0 - dt*beta
    return du
end
function d2step(u::Array{Float64,2}, s::Array{Float64,1})
    n = size(u)[2]
    ddu = zeros(3,9,n)
    for i = 1:n
    	x, y, z = u[1,i], u[2,i], u[3,i]
    	ddu[:,2,i] = [0., 0., -dt]
    	ddu[:,3,i] = [0., dt, 0.]
    	ddu[:,6,i] = [dt, 0., 0.]
    	ddu[:,8,i] = [-dt, 0., 0.]
    end
    return ddu
end
function d2step(u::Array{Float64,1}, s::Array{Float64,1})
    ddu = zeros(3,9)
    x, y, z = u[1], u[2], u[3]
    ddu[:,2] = [0., 0., -dt]
    ddu[:,3] = [0., dt, 0.]
    ddu[:,6] = [dt, 0., 0.]
    ddu[:,8] = [-dt, 0., 0.]
    return ddu
end
function pert(u::Array{Float64,2},s::Array{Float64,1})
    n = size(u)[2]
    x = view(u,1,:)
    pp = zeros(3, n)
    @. pp[2,:] = dt*x
    return pp
end
function pert(u::Array{Float64,1},s::Array{Float64,1})
    return [0, dt*u[1], 0]
end
function dpert(u::Array{Float64,1},s::Array{Float64,1})
    x, y, z = u[1], u[2], u[3]
    sigma, rho, beta = s
    dpp = zeros(3,3)
    dpp[2,:] .= dt
    deno  = 1/(dt^3*sigma*(beta*(rho - z - 1) + x*(y - x)) + dt^2*(sigma*(-rho + beta + z + 1) + beta + x^2) - dt*(sigma + beta + 1) + 1)
    du11 = (dt^2*beta + dt^2*x^2 - dt*beta - dt + 1)*deno
    du12 = (dt^2*sigma*beta - dt*sigma)*deno 
    du13 = (-dt^2*sigma*x)*deno
    dpp[2,1] *= du11
    dpp[2,2] *= du12
    dpp[2,3] *= du13
   
    return dpp
end
function dpert(u::Array{Float64,2},s::Array{Float64,1})
    x, y, z = u[1,:], u[2,:], u[3,:]
    sigma, rho, beta = s
    n = size(u)[2]
    dpp = zeros(3,3,n)
    dpp[2,:,:] .= dt
    
    deno  = @. 1/(dt^3*sigma*(beta*(rho - z - 1) + x*(y - x)) + dt^2*(sigma*(-rho + beta + z + 1) + beta + x^2) - dt*(sigma + beta + 1) + 1)
    du11 = @. (dt^2*beta + dt^2*x^2 - dt*beta - dt + 1)*deno
    du12 = @. (dt^2*sigma*beta - dt*sigma)*deno 
    du13 = @. (-dt^2*sigma*x)*deno
    dpp[2,1,:] .*= du11
    dpp[2,2,:] .*= du12
    dpp[2,3,:] .*= du13

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



