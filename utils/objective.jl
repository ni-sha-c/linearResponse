include("../examples/baker.jl")
using PyPlot
function obj_fun(x,y)
	return sin(x)*cos(y)	
end
function plot_obj_fun()

	n_xpts, n_ypts = 101, 101
	x_pts = LinRange(0.,2*pi, n_xpts)
	y_pts = LinRange(0.,2*pi, n_ypts)
	x_g = repeat(x_pts, 1, n_ypts) 
	y_g = repeat(y_pts', n_xpts, 1)
	J_g = obj_fun.(x_g, y_g)
	
	fig, ax = subplots(1,1)
	cm = plt.get_cmap("coolwarm")
	cplot = ax.contourf(x_g, y_g, J_g, 
				n_xpts, cmap=cm)
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	cbar = fig.colorbar(cplot, ax=ax)
	cbar.ax.tick_params(labelsize=28)
	ax.axis("scaled")

end
