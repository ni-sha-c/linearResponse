include("../examples/lorenz.jl")
using Test
#using PyPlot
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
function test_pert()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3)
	sp = s .+ [0., eps, 0.]
	sm = s .- [0., eps, 0.]
	
	pp_fd = (next(u, sp) - next(u, sm))/(2*eps)
	pp = pert(u, s)

	@show pp_fd, pp
	@test pp_fd ≈ pp atol=1.e-8
end
function test_pert_arr()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3)
	sp = s .+ [0., eps, 0.]
	sm = s .- [0., eps, 0.]
	
	pp_fd = (next(u,sp) - next(u, sm))/(2*eps)
	u = reshape(u, 3, 1)
	pp = pert(u, s)[:,1]

	@show pp_fd, pp
	@test pp_fd ≈ pp atol=1.e-8
end
function test_dpert_next()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3)
	u1px = u .+ eps*[1,0,0]
	u1mx = u .- eps*[1,0,0]
	u1py = u .+ eps*[0,1,0]
	u1my = u .- eps*[0,1,0]
	u1pz = u .+ eps*[0,0,1]
	u1mz = u .- eps*[0,0,1]

	dpert_dx_fd = (pert(u1px,s) - pert(u1mx,s))/(2*eps) 
	dpert_dy_fd = (pert(u1py,s) - pert(u1my,s))/(2*eps) 
	dpert_dz_fd = (pert(u1pz,s) - pert(u1mz,s))/(2*eps)
	dpp = dpert_next(u, s)

	@show dpp[:,1] .- dpert_dx_fd
	@show dpp[:,2] .- dpert_dy_fd
	@show dpp[:,3] .- dpert_dz_fd

	@test dpp[:,1] ≈ dpert_dx_fd atol=1.e-8
	@test dpp[:,2] ≈ dpert_dy_fd atol=1.e-8
	@test dpp[:,3] ≈ dpert_dz_fd atol=1.e-8
end
function test_dpert_next_arr()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3,1)
	u1px = u[:,1] .+ eps*[1,0,0]
	u1mx = u[:,1] .- eps*[1,0,0]
	u1py = u[:,1] .+ eps*[0,1,0]
	u1my = u[:,1] .- eps*[0,1,0]
	u1pz = u[:,1] .+ eps*[0,0,1]
	u1mz = u[:,1] .- eps*[0,0,1]

	dpert_dx_fd = (pert(u1px,s) - pert(u1mx,s))/(2*eps) 
	dpert_dy_fd = (pert(u1py,s) - pert(u1my,s))/(2*eps) 
	dpert_dz_fd = (pert(u1pz,s) - pert(u1mz,s))/(2*eps)
	dpp = dpert_next(u, s)

	@show dpp[:,1,1] .- dpert_dx_fd
	@show dpp[:,2,1] .- dpert_dy_fd
	@show dpp[:,3,1] .- dpert_dz_fd

	@test dpp[:,1,1] ≈ dpert_dx_fd atol=1.e-8
	@test dpp[:,2,1] ≈ dpert_dy_fd atol=1.e-8
	@test dpp[:,3,1] ≈ dpert_dz_fd atol=1.e-8
end
function test_dmag_flow()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3)
	dmag_dx = (mag_flow(u + eps*[1.0,0,0],s) - 
			   mag_flow(u - eps*[1.0,0,0],s))/(2*eps)
	dmag_dy = (mag_flow(u + eps*[0,1.0,0],s) - 
			   mag_flow(u - eps*[0,1.0,0],s))/(2*eps)
	dmag_dz = (mag_flow(u + eps*[0,0,1.0],s) - 
			   mag_flow(u - eps*[0,0,1.0],s))/(2*eps)

	dmag = dmag_flow(u, s)
	@test dmag[1] ≈ dmag_dx atol=1.e-8
	@test dmag[2] ≈ dmag_dy atol=1.e-8
	@test dmag[3] ≈ dmag_dz atol=1.e-8
end
function test_dmag_flow_arr()
	eps = 1.e-6
	s = [10., 28., 8.0/3]
	u = rand(3,1)
	dmag_dx = (mag_flow(u .+ eps*reshape([1.0,0,0],3,1),s) - 
			   mag_flow(u .- eps*reshape([1.0,0,0],3,1),s))/(2*eps)
	dmag_dy = (mag_flow(u .+ eps*reshape([0,1.0,0],3,1),s) - 
			   mag_flow(u .- eps*reshape([0,1.0,0],3,1),s))/(2*eps)
	dmag_dz = (mag_flow(u .+ eps*reshape([0,0,1.0],3,1),s) - 
			   mag_flow(u .- eps*reshape([0,0,1.0],3,1),s))/(2*eps)

	dmag = dmag_flow(u, s)
	@test dmag[1,:] ≈ dmag_dx atol=1.e-8
	@test dmag[2,:] ≈ dmag_dy atol=1.e-8
	@test dmag[3,:] ≈ dmag_dz atol=1.e-8
end
function test_dunit_flow()
	eps = 1.e-6
	s = [10., 28., 8/3]
	u = rand(3)
	upx = u .+ eps*[1.0,0,0]
	umx = u .- eps*[1.0,0,0]
	upy = u .+ eps*[0,1.0,0]
	umy = u .- eps*[0,1.0,0]
	upz = u .- eps*[0,0,1.0]
	umz = u .- eps*[0,0,1.0]
	dfdx = (flow(upx,s)/mag_flow(upx,s) .- 
			flow(umx,s)/mag_flow(umx,s))/(2*eps) 	
	dfdy = (flow(upy,s)/mag_flow(upy,s) .- 
			flow(umy,s)/mag_flow(umy,s))/(2*eps) 	
	dfdz = (flow(upz,s)/mag_flow(upz,s) .- 
			flow(umz,s)/mag_flow(umz,s))/(2*eps) 	
	df = dunit_flow(u, s)
	@show dfdx, df[:,1]

	@test df[:,1] ≈ dfdx atol=1.e-8
	@test df[:,2] ≈ dfdy atol=1.e-8
	@test df[:,3] ≈ dfdz atol=1.e-8

end
