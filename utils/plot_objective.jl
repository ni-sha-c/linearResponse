include("../examples/baker.jl")

include("get_objective.jl")
using PyPlot
using JLD
function plot_obj_erg_avg()
	X = load("../data/obj_erg_avg/cos4y_s4.jld")
	s4 = X["s4"]
	J = X["J"]
	fig, ax = subplots(1,1)
	ax.plot(s4, J, ".", ms=4.0)
	ax.grid(true)
	ax.set_xlim([minimum(s4), maximum(s4)])
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
end
function plot_obj_fun()

	n_xpts, n_ypts = 100, 100
	eps = 1.e-8
	x_pts = LinRange(0.,2*pi-eps, n_xpts)
	y_pts = LinRange(0.,2*pi-eps, n_ypts)
	x_g = repeat(x_pts, 1, n_ypts) 
	y_g = repeat(y_pts', n_xpts, 1)
	J_g = zeros(n_xpts, n_ypts)
	
	nSteps = 10000
	s = zeros(4)
	s[2] = 1.0
	for n = 1:nSteps
		J_g .+= obj_fun(x_g, y_g)/nSteps
		for j = 1:n_ypts, i = 1:n_xpts
			x_g[i,j], y_g[i,j] = step(
						[x_g[i,j], y_g[i,j]],
						s, 1)[:,end]
		end
	end
	@show minimum(J_g), maximum(J_g)
	fig, ax = subplots(1,1)
	cm = plt.get_cmap("coolwarm")
	J_g = J_g[:]
	clrs = cm(abs.(J_g)/maximum(abs.(J_g)))
	cplot = ax.scatter(x_g[:], y_g[:], 
					   c=clrs, cmap=cm)
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	#cbar = fig.colorbar(cplot, ax=ax)
	#cbar.ax.tick_params(labelsize=28)
	ax.axis("scaled")

end
