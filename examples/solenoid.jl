function cart_to_cyl(x,y)
    return [sqrt(x*x + y*y), (atan(y, x) + 2π)%(2π)]
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
	t1 = (2*t + 2π*s2*sin(2*t)) % (2π)
	z1 = z/s1 + sin(t)/2
	x1, y1 = cyl_to_cart(r1,t1)
    return [x1, y1, z1]
end
function dstep(u::Array{Float64,1},s::Array{Float64,1})
    du = zeros(2,2)
	s0, s1, s2 = s
    x, y, z = u
	r, t = cart_to_cyl(x,y)
	r1 = s0 + (r - s0)/s1 + cos(t)/2
	t1 = 2*t + 2π*s2*sin(2*t)
	z1 = z/s1 + sin(t)/2

	drt1drt = [1/s1 -sin(t)/2;0 2 + 4π*s2*cos(2*t)]
	dz1dt = cos(t)/2 
	dz1dz = 1/s1

	drtdxy = dcart_to_cyl(x,y)
	dtdx, dtdy = drtdxy[2,1], drtdxy[2,2]
	dxy1dxy = dcyl_to_cart(r1,t1)*drt1drt*drtdxy
	dz1 = [dz1dt*dtdx dz1dt*dtdy dz1dz]
	
	du = [dxy1dxy [0;0]]
	du = [du;dz1]
	return du
end
function pert(u::Array{Float64,1}, s::Array{Float64,1},
			  p::Int64)
	x, y, z = u
	s0, s1, s2 = s
 	r, t = cart_to_cyl(x,y)	 
    r1 = s0 + (r - s0)/s1 + cos(t)/2
	t1 = 2*t + 2π*s2*sin(2*t)
	z1 = z/s1 + sin(t)/2
	x1, y1 = cyl_to_cart(r1,t1)
	dxyz1drtz1 = [dcyl_to_cart(r1,t1) [0;0]]
	dxyz1drtz1 = [dxyz1drtz1; [0 0 1]] 

	if p==1
		drtz1ds = [1-1/s1; 0; 0]
    elseif p==2
		drtz1ds = [-(r - s0)/s1^2; 0; -z/s1^2]
    elseif p==3
		drtz1ds = [0; 2π*sin(2*t); 0]  
    else
    	println("parameter perturbation is only defined 
    			for indices 1 through 3")
    	return 0
    end
	return dxyz1drtz1*drtz1ds
end
