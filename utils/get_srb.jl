include("../examples/baker.jl")
using JLD
using SharedArrays
using Distributed
function smooth_indicator(x1,x2,a1,a2,h1,h2)
	f1, f2 = 0., 0.
	if (-h1/2 <= x1 - a1 <= h1/2)
		f1 = cos(pi*(x1-a1)/h1)
	end
	if	(-h1/2 <= (x1 - 2*pi) - a1 <= h1/2) 
		f1 = cos(pi*(x1-a1-2*pi)/h1)
	end
	if (-h1/2 <= (x1 + 2*pi) - a1 <= h1/2) 
		f1 = cos(pi*(x1-a1+2*pi)/h1)
	end
	if (-h2/2 <= x2 - a2 <= h2/2) 
		f2 = cos(pi*(x2-a2)/h2)
	end
	if	(-h2/2 <= (x2 - 2*pi) - a2 <= h2/2) 
		f2 = cos(pi*(x2-a2-2*pi)/h2)
	end
	if	(-h2/2 <= (x2 + 2*pi) - a2 <= h2/2)
		f2 = cos(pi*(x2-a2+2*pi)/h2)
	end
	return f1*f2
end
function compute_indicator_density(s)
	n_xbins, n_ybins = 200, 200
	dx, dy = 2*pi/n_xbins, 2*pi/n_ybins

	rho = zeros(n_xbins, n_ybins)
	rho .= 0.
	n_step = 10000
	n_spl = nprocs() - 1 
	n_rep = 10000
	rho_proc = SharedArray{Float64}(n_xbins*n_ybins,
					 n_spl)

	for i = 1:n_rep
		rho_proc .= 0.
		t = @distributed for i = 1:n_spl
			u = 2*pi*rand(2)
			u_trj = step(u, s, n_step-1)
			x, y = view(u_trj,1,:),view(u_trj,2,:)
			x_ind = floor.(Int64, x/dx) .+ 1
			y_ind = floor.(Int64, y/dy) .+ 1
			locs = (y_ind .- 1)*n_xbins .+ 
				x_ind
			for l in locs
				rho_proc[l,i] += 
				1.0/n_step/n_spl/n_rep
			end
		end
		wait(t)
		for k = 1:n_spl
			rho[:] .+= rho_proc[:,k]/dx/dy
		end
		
	end
	save(string("../data/SRB_dist/ind_dist_", 
				"$s","_.jld"), "rho", rho)
end
function get_dist()
	s = zeros(4)
	s[3] = 1.0
	compute_indicator_density(s)
end
