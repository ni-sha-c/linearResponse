function cart_to_cyl(x,y)
    return [sqrt(x*x + y*y), (atan2(y, x) + 2π)%(2π)]
end
function cyl_to_cart(r,t)
    return [r*cos(t), r*sin(t)] 
end
function dcart_to_cyl(x,y)
	r2 = x*x + y*y 
	r = sqrt(r2)
	return [x/r y/r; -y/r2 x/r2]
end
function dcyl_to_cart(r,t)
	ct, st = cos(t), sin(t)
	return [ct -r*st;st r*ct]	
end
function step(x, s, n)
    x_trj = zeros(3, n+1)
    x_trj[:,1] = x
	s0, s1, s2 = s[1], s[2], s[3]
    for i = 2:n+1
   		x_trj[:,i] = next(x_trj[:,i-1],s)	
	end
    return x_trj
end
function next(u, s)
    x, y, z = u
	s0, s1, s2 = s
	r, t = cart_to_cyl(x, y)
	r1 = s0 + (r - s0)/s1 + cos(t)/2
	t1 = 2*t + 2π*s2*sin(2*t)
	z1 = z/s1 + sin(t)/2
	x1, y1 = cyl_to_cart(r1,t1)
    return [x1, y1, z1]
end

function dstep(u::Array{Float64,1},s::Array{Float64,1})
    du = zeros(2,2)
	s0, s1, s2 = s
    x, y, z = u
	r, t = cart_to_cyl(x,y)
	dr1dr = 1/s1
	dr1dt = -sin(t)/2
	dt1dt = 2 + 4π*s2*cos(2*t)
	dz1dt = cos(t)/2 
	dz1dz = 1/s1

	return du
end

function dstep(u::Array{Float64,2},s::Array{Float64,1})
    n = size(u)[2]
    du = zeros(2,2,n)
    x, y = view(u,1,:),view(u,2,:)
    sx, sy = sin.(x), sin.(2*y)/2
    dsx, dsy = cos.(x), cos.(2*y)
    du[1,1,:] = 2 .+ s[1]*dsx .+ 
    			s[2]*sy.*dsx 
    du[1,2,:] = s[2]*sx.*dsy
    du[2,1,:] = s[3]*dsx.*sy
    du[2,2,:] = 0.5 .+ s[4]*dsy .+ 
    			s[3]*sx.*dsy
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
