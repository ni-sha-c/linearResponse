function step(x, s, n)
	x_trj = zeros(2, n+1)
	x_trj[:,1] = x
	for i = 2:n+1
		x, y = x_trj[1,i-1],x_trj[2,i-1]
		sx, sy = sin(x), sin(2*y)/2
		x1 = (2*x + 
			(s[1] + s[2]*sy)*sx) 
		y1 = (0.5*y + 
			(s[4] + s[3]*sx)*sy)   
		
		x_trj[1,i] = x < pi ? x1 : x1 - 2*pi 
		x_trj[2,i] = x < pi ? y1 : y1 + pi
		
		x_trj[1,i] = x_trj[1,i] % (2*pi)
		x_trj[2,i] = x_trj[2,i] % (2*pi)

	end
	return x_trj
end
function next(u, s)
	x, y = u[1], u[2]
	sx, sy = sin(x), sin(2*y)/2
	x1 = (2*x + 
			(s[1] + s[2]*sy)*sx) 
	y1 = (0.5*y + 
			(s[4] + s[3]*sx)*sy)   
		
	x_next = x < pi ? x1 : x1 - 2*pi 
	y_next = x < pi ? y1 : y1 + pi
		
	x_next = x_next % (2*pi)
	y_next = y_next % (2*pi)

	return [x_next, y_next]
end

function dstep(u::Array{Float64,1},s::Array{Float64,1})
	du = zeros(2,2)
	x, y = u[1], u[2]
	sx, sy = sin(x), sin(2*y)/2
	dsx, dsy = cos(x), cos(2*y)
	du[1,1] = 2 + s[1]*dsx + s[2]*sy*dsx 
	du[1,2] = s[2]*sx*dsy
	du[2,1] = s[3]*dsx*sy
	du[2,2] = 0.5 + s[4]*dsy + s[3]*sx*dsy
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
function pert(u, p)
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
