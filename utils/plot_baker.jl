include("../examples/baker.jl")
using PyPlot
function plot_one_step()
	n_xq, n_yq = 20, 40
	eps = 1.e-14
	x_q = LinRange(0, pi - eps, n_xq)	
	y_q = LinRange(0, 2*pi, n_yq)	
	#x_q = x_q[1:end-1]
	#y_q = y_q[1:end-1]
	#n_xq = n_xq - 1
	#n_yq = n_yq - 1
	cm = plt.get_cmap("winter")
	x_q = repeat(x_q, 1, n_yq)
	y_q = repeat(y_q', n_xq, 1)
	fig, ax = subplots(1,1)
	fig1, ax1 = subplots(1,1)
	
	cplot = ax.contour(x_q, y_q, y_q, n_yq, 
					   cmap=cm)
	cbar = fig.colorbar(cplot, ax=ax)
	x1_q = zeros(n_xq, n_yq)
	y1_q = zeros(n_xq, n_yq)

	for j = 1:n_yq, i = 1:n_xq
		x1_q[i,j], y1_q[i,j] = 
			step([x_q[i,j], y_q[i,j]]
				, zeros(4), 1)[:,end]
	end
	cplot1 = ax1.contour(x1_q, y1_q, y_q, n_yq, 
					   cmap=cm)
	x_q = x_q'
	y_q = y_q'
	x1_q = x1_q'
	y1_q = y1_q'
	ax.contour(x_q, y_q, x_q, n_xq, 
					   cmap=cm, vmin=0.,vmax=2*pi)
	ax1.contour(x1_q, y1_q, x_q, n_xq, 
					   cmap=cm, vmin=0., vmax=2*pi)




	#n_xq = n_xq + 1
	#n_yq = n_yq + 1
	x_q = LinRange(pi + eps, 2*pi, n_xq)	
	#x_q = x_q[1:end-1]
	y_q = LinRange(0., 2*pi, n_yq)
	#y_q = y_q[1:end-1]
	#n_xq = n_xq - 1
	#n_yq = n_yq - 1
	x_q = repeat(x_q, 1, n_yq)
	y_q = repeat(y_q', n_xq, 1)
	
	ax.contour(x_q, y_q, y_q, n_yq, 
				cmap=cm)
	x1_q = zeros(n_xq, n_yq)
	y1_q = zeros(n_xq, n_yq)

	for j = 1:n_yq, i = 1:n_xq
		x1_q[i,j], y1_q[i,j] = 
			step([x_q[i,j], y_q[i,j]]
				, zeros(4), 1)[:,end]
	end
	cplot1 = ax1.contour(x1_q, y1_q, y_q, n_yq, 
					   cmap=cm)
	cbar1 = fig1.colorbar(cplot, ax=ax1)
	x_q = x_q'
	y_q = y_q'
	x1_q = x1_q'
	y1_q = y1_q'
	ax.contour(x_q, y_q, x_q, n_xq, 
					   cmap=cm,vmin=0.,vmax=2*pi)
	ax1.contour(x1_q, y1_q, x_q, n_xq, 
					   cmap=cm, vmin=0., vmax=2*pi)

	ax1.xaxis.set_tick_params(labelsize=28)
	ax1.yaxis.set_tick_params(labelsize=28)

	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	cbar.ax.tick_params(labelsize=28)
	cbar1.ax.tick_params(labelsize=28)
	ax.axis("scaled")
	ax1.axis("scaled")
	ax.set_xlim([0.,2*pi])
	ax.set_ylim([0.,2*pi])
	ax1.set_xlim([0.,2*pi])
	ax1.set_ylim([0.,2*pi])
end
function plot_mapping_vs_params()
	n_pts = 101
	x_pts = LinRange(0, 2*pi, n_pts)	
	x_pts = x_pts[1:end-1]
	n_pts = n_pts - 1

	cm = plt.get_cmap("coolwarm")
	clrs = cm(x_pts/(2*pi))
	x_gr = repeat(x_pts, 1, n_pts)
	y_gr = x_gr'
	
	for n = 1:4
		s = zeros(4)
		s[n] = 1.0
	
		x1_gr = zeros(n_pts, n_pts)
		y1_gr = zeros(n_pts, n_pts)

		for j = 1:n_pts, i = 1:n_pts
			x1_gr[i,j], y1_gr[i,j] = 
				step([x_gr[i,j], y_gr[i,j]]
					, s, 1)[:,end]
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
		cbar.ax.tick_params(labelsize=28)
		ax1.axis("scaled")
		ax1.set_xlim([0.,2*pi])
		ax1.set_ylim([0.,2*pi])
	end
end












