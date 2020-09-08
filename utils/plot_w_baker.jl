include("../examples/baker.jl")
using LinearAlgebra
using PyPlot
using PyCall
function compute_w()
	s = [0.3, 0., 0.3, 0.0]
	n = 40
	m = 40
	x = LinRange(0.,2*pi,n)
	u = [[x[i] + 2*pi*rand(),x[j] + 2*pi*rand()] 
		 for i=1:n,j = 1:n][:]
	q = [rand(2) for i=1:n, j=1:n][:]
	w = [zeros(2) for i=1:n, j=1:n][:]
	D2 = [zeros(2,2) for i=1:n, j=1:n][:]
	n = length(u)
	z = zeros(n)
	pts = zeros(2, n)
	vecs = zeros(2, n)
	eps = 1.e-1
	fig, ax = subplots(1,1)
	for i = 1:m
		@show i
		u .= next.(u, Ref(s))
		q .= pushforward.(u, q, Ref(s))
		z .= norm.(q)
		q .= q./z
	end
	for i = 1:m
		@show i
		D2 .= pushforward_second_order.(u, q, Ref(s))	
		w .= pushforward.(u, w, Ref(s))
		w .+= tensordot(D2, q)
		q .= pushforward.(u, q, Ref(s))
		z .= norm.(q)
		q .= q./z
		w ./= (z.*z)
		w .= w .- dot.(w, q).*q
		u .= next.(u, Ref(s))
		if rem(i, 10) == 1
			pts .= hcat(u...)
			vecs .= hcat(w...)
			x_pts = [pts[1,:] - eps*vecs[1,:] pts[1,:] + eps*vecs[1,:]]'

			y_pts =	[pts[2,:] - eps*vecs[2,:] pts[2,:] + eps*vecs[2,:]]' 
			ax.clear()
			ax.set_xlim([0,2*pi])
			ax.set_ylim([0,2*pi])
			ax.set_xlabel(L"$x_1$", fontsize=30)
			ax.set_ylabel(L"$x_2$", fontsize=30)
			ax.xaxis.set_tick_params(labelsize=30)
			ax.yaxis.set_tick_params(labelsize=30)
			ax.axis("scaled")

			ax.plot(x_pts, y_pts, "r")
			savefig(string("plots/w_n_",i,".png"))
		end
	end 
	return pts, vecs

end
function tensordot(A, b)
	return [Ai*b[i] for (i, Ai) in enumerate(A)]
end
function pushforward_second_order(u, v1, s)
	d2u = d2step(u, s)
	d2u_v1 = reshape([dot(d2u[:,i], v1) for i=1:4], 2, 2)
	return d2u_v1
end

function pushforward(u, q, s)
	du = dstep(u, s)
	q = du*q
	return q
end
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
	n_xpts, n_ypts = 75, 75
	x_pts = LinRange(0, 2*pi, n_xpts)	
	#x_pts = x_pts[1:end-1]
	#n_xpts = n_xpts - 1

	y_pts = LinRange(0, 2*pi, n_ypts)	
	#y_pts = y_pts[1:end-1]
	#n_ypts = n_ypts - 1


	cm = plt.get_cmap("winter")
	x_gr = repeat(x_pts, 1, n_ypts)
	y_gr = repeat(y_pts', n_xpts, 1)
	x_gr = x_gr[:]
	y_gr = y_gr[:]
	n_gr = length(x_gr)

	#n_xpts = n_xpts + 1
	#n_ypts = n_ypts + 1
	eps = 1.e-1
	x_pts = LinRange(eps, 2*pi, n_xpts)	
	#x_pts = x_pts[1:end-1]
	#n_xpts = n_xpts - 1

	y_pts = LinRange(eps, 2*pi, n_ypts)	
	#y_pts = y_pts[1:end-1]
	#n_ypts = n_ypts - 1
	
	x_gr1 = repeat(x_pts, 1, n_ypts)
	y_gr1 = repeat(y_pts', n_xpts, 1)
	
	x_gr1 = x_gr1[:]
	y_gr1 = y_gr1[:]


	for n = 1:4
		s = zeros(4)
		s[n] = 1.0
	
		x1_gr = zeros(n_gr)
		y1_gr = zeros(n_gr)

		for i = 1:n_gr
			x1_gr[i], y1_gr[i] = 
				step([x_gr[i], y_gr[i]]
					 , s, 1)[:,end]
		end
		# Horizontal lines
		fig1, ax1 = subplots(1,1)
		clrs = cm(y_gr/(2*pi))
		sp1 = ax1.scatter(x1_gr, y1_gr, s=40.0,
						  c=clrs, marker="_", 
					   	  cmap=cm)
		#cbar = fig.colorbar(sp1, ax=ax1)

	
		# Vertical lines
		x1_gr1 = zeros(n_gr)
		y1_gr1 = zeros(n_gr)

		for i = 1:n_gr
			x1_gr1[i], y1_gr1[i] = 
				step([x1_gr1[i], y1_gr1[i]]
					 , s, 1)[:,end]
		end


		clrs = cm(x_gr1/(2*pi))
		ax1.scatter(x_gr1, y_gr1, s=40.0, 
					marker="|",c=clrs, cmap=cm)
		ax1.xaxis.set_tick_params(labelsize=28)
		ax1.yaxis.set_tick_params(labelsize=28)
		ax1.axis("scaled")
		ax1.set_xlim([0.,2*pi])
		ax1.set_ylim([0.,2*pi])
		ax1.set_title("\$s_$n = 1\$", fontsize=28)
		
	end
end












