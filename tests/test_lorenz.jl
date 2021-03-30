include("../examples/lorenz.jl")
using Test
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
function test_dstep()
	eps = 1.e-8
	s = [10., 28., 8.0/3]
	u = rand(3)
	u1px = u .+ eps*[1,0,0]
	u1mx = u .- eps*[1,0,0]
	u1py = u .+ eps*[0,1,0]
	u1my = u .- eps*[0,1,0]
	u1pz = u .+ eps*[0,0,1]
	u1mz = u .- eps*[0,0,1]

	du_dx_fd = (next(u1px,s) - next(u1mx,s))/(2*eps) 
	du_dy_fd = (next(u1py,s) - next(u1my,s))/(2*eps) 
	du_dz_fd = (next(u1pz,s) - next(u1mz,s))/(2*eps)

	du = dstep(u, s)
	du_dx = du[:,1]
	du_dy = du[:,2]
	du_dz = du[:,3]

	@show abs.(du_dx_fd .- du_dx)
	@show abs.(du_dy_fd .- du_dy)
	@show abs.(du_dz_fd .- du_dz)
	@test du_dx_fd ≈ du_dx atol=1.e-8
	@test du_dy_fd ≈ du_dy atol=1.e-8
	@test du_dz_fd ≈ du_dz atol=1.e-8
end
function test_dstep_arr()
	eps = 1.e-8
	s = [10., 28., 8.0/3]
	u = rand(3,1)
	u1px = u[:,1] .+ eps*[1,0,0]
	u1mx = u[:,1] .- eps*[1,0,0]
	u1py = u[:,1] .+ eps*[0,1,0]
	u1my = u[:,1] .- eps*[0,1,0]
	u1pz = u[:,1] .+ eps*[0,0,1]
	u1mz = u[:,1] .- eps*[0,0,1]

	du_dx_fd = (next(u1px,s) - next(u1mx,s))/(2*eps) 
	du_dy_fd = (next(u1py,s) - next(u1my,s))/(2*eps) 
	du_dz_fd = (next(u1pz,s) - next(u1mz,s))/(2*eps)

	du = dstep(u, s)
	@show du
	du_dx = du[:,1,1]
	du_dy = du[:,2,1]
	du_dz = du[:,3,1]

	@test du_dx_fd ≈ du_dx atol=1.e-8
	@test du_dy_fd ≈ du_dy atol=1.e-8
	@test du_dz_fd ≈ du_dz atol=1.e-8
end
function test_d2step()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3)
	u1px = u .+ eps*[1,0,0]
	u1mx = u .- eps*[1,0,0]
	u1py = u .+ eps*[0,1,0]
	u1my = u .- eps*[0,1,0]
	u1pz = u .+ eps*[0,0,1]
	u1mz = u .- eps*[0,0,1]

	ddu_dx_fd = (dstep(u1px,s) - dstep(u1mx,s))/(2*eps) 
	ddu_dy_fd = (dstep(u1py,s) - dstep(u1my,s))/(2*eps) 
	ddu_dz_fd = (dstep(u1pz,s) - dstep(u1mz,s))/(2*eps)

	ddu = d2step(u, s)
	@show ddu[:,:,1] .- ddu_dx_fd
	@show ddu[:,:,2] .- ddu_dy_fd
	@show ddu[:,:,3] .- ddu_dz_fd

	@test ddu[:,:,1] ≈ ddu_dx_fd atol=1.e-8
	@test ddu[:,:,2] ≈ ddu_dy_fd atol=1.e-8
	@test ddu[:,:,3] ≈ ddu_dz_fd atol=1.e-8

end
function test_d2step_arr()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3,1)
	u1px = u[:,1] .+ eps*[1,0,0]
	u1mx = u[:,1] .- eps*[1,0,0]
	u1py = u[:,1] .+ eps*[0,1,0]
	u1my = u[:,1] .- eps*[0,1,0]
	u1pz = u[:,1] .+ eps*[0,0,1]
	u1mz = u[:,1] .- eps*[0,0,1]

	ddu_dx_fd = (dstep(u1px,s) - dstep(u1mx,s))/(2*eps) 
	ddu_dy_fd = (dstep(u1py,s) - dstep(u1my,s))/(2*eps) 
	ddu_dz_fd = (dstep(u1pz,s) - dstep(u1mz,s))/(2*eps)

	ddu = d2step(u, s)
	@show ddu[:,:,1,1] .- ddu_dx_fd
	@show ddu[:,:,2,1] .- ddu_dy_fd
	@show ddu[:,:,3,1] .- ddu_dz_fd

	@test ddu[:,:,1,1] ≈ ddu_dx_fd atol=1.e-8
	@test ddu[:,:,2,1] ≈ ddu_dy_fd atol=1.e-8
	@test ddu[:,:,3,1] ≈ ddu_dz_fd atol=1.e-8

end
