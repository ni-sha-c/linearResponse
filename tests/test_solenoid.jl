include("../examples/solenoid.jl")
using Test
using Zygote
#using PyPlot
function test_attractor()
	s = [1.0, 4.0, 0.1]
	u = rand(3)
	n_runup = 2000
	for i = 1:n_runup
		u = next(u,s)
	end
	@show u
	n = 50000
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
	eps = 1.e-6
	s = [1., 4., 0.]
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
function test_d2step()
# d2step has to be defined in the main script as below, as opposed to writing a function in solenoid.jl, because we don't want to do source-to-source AD every time the second derivative is needed.
	eps = 1.e-6
	s = [1.0, 4.0,0.1]
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

	ddu(x) = jacobian(x -> dstep(x,s), x)
	ddu_ad = ddu(u)[1]
	ddu_dx_ad = reshape(ddu_ad[:,1],3,3)
	ddu_dy_ad = reshape(ddu_ad[:,2],3,3)
	ddu_dz_ad = reshape(ddu_ad[:,3],3,3)

	@show ddu_dx_ad .- ddu_dx_fd
	@show ddu_dy_ad .- ddu_dy_fd
	@show ddu_dz_ad .- ddu_dz_fd

	@test ddu_dx_ad ≈ ddu_dx_fd atol=1.e-8
	@test ddu_dy_ad ≈ ddu_dy_fd atol=1.e-8
	@test ddu_dz_ad ≈ ddu_dz_fd atol=1.e-8

end
function test_pert()
	eps = 1.e-6
	s = [1., 4., 0.]
	u = rand(3)
	for i = 1:3
		eps_arr = zeros(3)
		eps_arr[i] = eps 
		sp = s .+ eps_arr
		sm = s .- eps_arr
	
		pp_fd = (next(u, sp) - next(u, sm))/(2*eps)
		pp = pert(u, s, i)
		@show i
		@show pp_fd, pp
		@test pp_fd ≈ pp atol=1.e-8
	end
end

