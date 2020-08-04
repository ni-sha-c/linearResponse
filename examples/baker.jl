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
		
		x_trj[1,i] = x < pi ? x1 : 4*pi - x1
		x_trj[2,i] = x < pi ? y1 : y1 + pi

	end
	return x_trj
end
function dstep(u,s)
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
