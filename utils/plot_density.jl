using PyPlot
using JLD
include("get_srb.jl")
function plot_smooth_indicator()
	fig = figure() 
	ax = fig.add_subplot(111, projection="3d")
	eps = 1.e-8
	n_xpts, n_ypts = 100,100
	x1 = LinRange(0., 2*pi-eps, n_xpts)
	x2 = LinRange(0., 2*pi-eps, n_ypts)
	h1, h2 = 1.0, 1.0
	x1_gr = repeat(x1, 1, n_ypts)
	x2_gr = repeat(x2', n_xpts, 1)
	n_ctrs = 10
	a = 2*pi*rand(n_ctrs)
	b = 2*pi*rand(n_ctrs)
	x3_gr = zeros(n_xpts, n_ypts, n_ctrs)
	for (i, a1) = enumerate(a)
		x3_gr[:,:,i] = smooth_indicator.(x1_gr,
						x2_gr, a1, b[i], h1, h2)
		ax.plot_surface(x1_gr, x2_gr, x3_gr[:,:,i], 
				label="\$ a = $a1 \$",cmap="Blues",
				vmin=0.,vmax=1.0)
	end
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax.zaxis.set_tick_params(labelsize=28)
	ax.set_xlabel("x",fontsize=28)
	ax.set_ylabel("y",fontsize=28)
	ax.set_zlabel("z",fontsize=28)
	ax.set_title("\$ a = $(a[1]) \$",fontsize=28)
	#lgnd = fig.legend(fontsize=28)
end
function plot_density()
	#assumes get_dist() has been run
	s = zeros(4)
	s[1] = 1.0
	X = load(string("../data/SRB_dist/ind_dist_", 
				   "$s","_.jld"))
	rho = X["rho"]
	n_x, n_y = size(rho)
	x = LinRange(0.,2*pi,n_x)
	y = LinRange(0.,2*pi,n_y)
	x_g = repeat(x, 1, n_y)
	y_g = repeat(y', n_x, 1)

	fig, ax = subplots(1,1)
	cplot = ax.contour(x_g, y_g, rho, n_x,
					   cmap="Blues",vmin=0.,
					   vmax=1.0)
	cbar = fig.colorbar(cplot,ax=ax)
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	cbar.ax.tick_params(labelsize=28)
	ax.set_title("\$ s_1 = $(s[1]) \$",fontsize=28)
end
