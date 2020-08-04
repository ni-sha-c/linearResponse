include("../examples/baker.jl")
using PyPlot
function plot_one_step()
	n_pts = 50
	x_pts = LinRange(0, 2*pi, n_pts)	
	#x_pts = x_pts[1:end-1]
	#n_pts = n_pts - 1
	cm = plt.get_cmap("coolwarm")
	clrs = cm(x_pts/(2*pi))
	x_gr = repeat(x_pts, 1, n_pts)
	y_gr = x_gr'

	# Horizontal lines
	fig, ax = subplots(1,1)
	cplot = ax.contour(x_gr, y_gr, y_gr, n_pts, 
					   cmap=cm)
	cbar = fig.colorbar(cplot, ax=ax)

	# Vertical lines
	ax.contour(y_gr, x_gr, y_gr, n_pts, cmap=cm)
	
	# Iterate once
	x1_gr = zeros(n_pts, n_pts)
	y1_gr = zeros(n_pts, n_pts)

	for j = 1:n_pts, i = 1:n_pts
		x1_gr[i,j], y1_gr[i,j] = 
			step([x_gr[i,j], y_gr[i,j]]
				, zeros(4), 1)[:,end]
	end
	# Horizontal lines
	fig1, ax1 = subplots(1,1)
	cplot1 = ax1.contour(x1_gr, y1_gr, y_gr, n_pts, 
					   cmap=cm)
	cbar = fig.colorbar(cplot1, ax=ax1)

	# Vertical lines
	ax1.contour(x1_gr', y1_gr', y_gr, n_pts, cmap=cm)
	ax1.xaxis.set_tick_params(labelsize=28)
	ax1.yaxis.set_tick_params(labelsize=28)

	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	cbar.ax.tick_params(labelsize=28)
	ax.axis("scaled")
	ax1.axis("scaled")
	ax.set_xlim([0.,2*pi])
	ax.set_ylim([0.,2*pi])
	ax1.set_xlim([0.,2*pi])
	ax1.set_ylim([0.,2*pi])
end
function plot_mapping_vs_params()













end
