include("../examples/lorenz.jl")
using PyPlot
function test_attractor()
	s = [10., 28., 8.0/3]
	u = rand(3)
	n_runup = 20000
	for i = 1:n_runup
		u = next(u,s)
	end
	@show u
	n = 500000
	u_trj = step(u, s, n)
	u_trj = u_trj[:,1:10:end]	
	x, y, z = view(u_trj,1,:), view(u_trj,2,:), view(u_trj,3,:)
	fig, ax = subplots(1,1)
	ax.plot(x,y,".",ms=2)
	ax.set_xlabel("x",fontsize=36)
	ax.set_ylabel("y",fontsize=36)
	ax.tick_params(axis="both",labelsize=36)
	fig, ax = subplots(1,1)
	ax.plot(x,z,".",ms=2)
	ax.set_xlabel("x",fontsize=36)
	ax.set_ylabel("z",fontsize=36)
	ax.tick_params(axis="both",labelsize=36)	

end
