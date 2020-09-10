include("../examples/baker.jl")
using LinearAlgebra
using PyPlot
using PyCall
@pyimport matplotlib.collections as coll 
@pyimport matplotlib.colors as cs 
@pyimport matplotlib.cm as cm 
function plot_y(m)
	s = [0.3, 0., 0.3, 0.0]
	n = 400
	x = LinRange(0.,2*pi,n)
	u = [[x[i] + 2*pi*rand(),x[j] + 2*pi*rand()] 
		 for i=1:n,j = 1:n][:]
	q = [rand(2) for i=1:n, j=1:n][:]
	w = [zeros(2) for i=1:n, j=1:n][:]
	y = [zeros(2) for i=1:n, j=1:n][:]
	v = [zeros(2) for i=1:n, j=1:n][:]
	X = [zeros(2) for i=1:n, j=1:n][:]

	D2 = [zeros(2,2) for i=1:n, j=1:n][:]
	DX = [zeros(2,2) for i=1:n, j=1:n][:]
	n = length(u)
	z = zeros(n)
	a = zeros(n)
	delta = zeros(n)
	ny = zeros(n)
	segments = zeros(n, 2, 2)
	eps = 1.e-1
	filepath = string("/home/nishac/Research/PhDThesis/papers/",
					  "Decomposition-of-Linear-Response/figs/")
	for i = 1:20
		@show i
		u .= next.(u, Ref(s))
		q .= pushforward.(u, q, Ref(s))
		z .= norm.(q)
		q .= q./z
	end


	for i = 1:m
		@show i
		D2 .= pushforward_second_order.(u, q, Ref(s))	
		y .= pushforward.(u, y, Ref(s))
		y .+= tensordot(D2, v)
		X .= pert.(u, 1) .+ pert.(u, 3) 
		DX .= dpert.(u,1) .+ dpert.(u,3)
		v .= pushforward.(u, v, Ref(s))
		v .+= X
		w .= pushforward.(u, w, Ref(s))
		w .+= tensordot(D2, q)

		q .= pushforward.(u, q, Ref(s))
		z .= norm.(q)
		q .= q./z
	
		a .= dot.(v, q)
		v .-= a.*q

		w ./= (z.*z)
		w .= w .- dot.(w, q).*q
		y ./= z	
		y .+= -a.*w .+ tensordot(DX, q)
		delta .= dot.(y, q) .+ dot.(v, w)
		y .-= delta.*q
		u .= next.(u, Ref(s))
	end
	qperp = [[-qi[2], qi[1]] for qi in q]
	ny = dot.(y, qperp) 
	segments .= create_line_colls(u, q, eps)
	lc = coll.LineCollection(segments, 
						cmap=plt.get_cmap("coolwarm"),
						norm=cs.Normalize(minimum(ny),
						maximum(ny)))
	lc.set_array(ny)
	lc.set_linewidth(2)
	
	fig, ax = subplots(1,1)
	ax.set_xlim([0,2*pi])
	ax.set_ylim([0,2*pi])
	ax.set_xlabel(L"$x_1$", fontsize=30)
	ax.set_ylabel(L"$x_2$", fontsize=30)
	ax.xaxis.set_tick_params(labelsize=30)
	ax.yaxis.set_tick_params(labelsize=30)

	ax.add_collection(lc)
	ax.axis("scaled")

	cbar = fig.colorbar(cm.ScalarMappable(
						norm=cs.Normalize(minimum(ny),
						maximum(ny)), 
					   cmap=plt.get_cmap("coolwarm")), ax=ax,
						orientation="horizontal",shrink=0.4,
						pad=0.1)

	cbar.ax.tick_params(labelsize=30)
	cbar.ax.xaxis.get_offset_text().set_fontsize(30)
	plt.tight_layout()
	#savefig(string(filepath,"w_n_",m,".png"))
	fig, ax = subplots(1,1)
	ax.set_xlim([0,2*pi])
	ax.set_ylim([0,2*pi])
	ax.set_xlabel(L"$x_1$", fontsize=30)
	ax.set_ylabel(L"$x_2$", fontsize=30)
	ax.xaxis.set_tick_params(labelsize=30)
	ax.yaxis.set_tick_params(labelsize=30)
	pts = hcat(u...)
	vecs = hcat(y...)
	x_pts = [pts[1,:] - eps*vecs[1,:] pts[1,:] + eps*vecs[1,:]]'
	y_pts = [pts[2,:] - eps*vecs[2,:] pts[2,:] + eps*vecs[2,:]]'

	ax.plot(x_pts, y_pts, "darkcyan")
	ax.axis("scaled")


	return segments, ny
end

function plot_v(m)
	s = [0.1, 0.1, 0., 0.1]
	n = 400
	x = LinRange(0.,2*pi,n)
	u = [[x[i] + 2*pi*rand(),x[j] + 2*pi*rand()] 
		 for i=1:n,j = 1:n][:]
	q = [rand(2) for i=1:n, j=1:n][:]
	v = [zeros(2) for i=1:n, j=1:n][:]
	X = [zeros(2) for i=1:n, j=1:n][:]
	n = length(u)
	z = zeros(n)
	a = zeros(n)
	segments = zeros(n, 2, 2)
	eps = 1.e-1
	eps1 = 1.e-4
	filepath = string("/home/nishac/Research/PhDThesis/papers/",
					  "Decomposition-of-Linear-Response/figs/")
	y_pos = 1.0
	for i = 1:20
		@show i
		u .= next.(u, Ref(s))
		q .= pushforward.(u, q, Ref(s))
		z .= norm.(q)
		q .= q./z
	end
	for i = 1:m
		@show i
		X .= pert.(u, 1) .+ pert.(u, 3) 
		v .= pushforward.(u, v, Ref(s))
		v .+= X
		q .= pushforward.(u, q, Ref(s))
		z .= norm.(q)
		q .= q./z
		a .= dot.(v, q)
		v .-= a.*q
		u .= next.(u, Ref(s))

	end
	inds = [abs(uk[2] - y_pos) < eps1 for uk in u]
	u, v, q = u[inds], v[inds], q[inds]
	#=
	a .= norm.(v)
	segments .= create_line_colls(u, q, eps)
	lc = coll.LineCollection(segments, 
						cmap=plt.get_cmap("coolwarm"),
						norm=cs.Normalize(minimum(a),
						maximum(a)))
	lc.set_array(a)
	lc.set_linewidth(2)
	
	fig, ax = subplots(1,1)
	ax.set_xlim([0,2*pi])
	ax.set_ylim([0,2*pi])
	ax.set_xlabel(L"$x_1$", fontsize=30)
	ax.set_ylabel(L"$x_2$", fontsize=30)
	ax.xaxis.set_tick_params(labelsize=30)
	ax.yaxis.set_tick_params(labelsize=30)

	ax.add_collection(lc)
	ax.axis("scaled")

	cbar = fig.colorbar(cm.ScalarMappable(
						norm=cs.Normalize(minimum(a),
						maximum(a)), 
					   cmap=plt.get_cmap("coolwarm")), ax=ax,
						orientation="horizontal",shrink=0.4,
						pad=0.1)

	cbar.ax.tick_params(labelsize=30)
	cbar.ax.xaxis.get_offset_text().set_fontsize(30)
	plt.tight_layout()
	#savefig(string(filepath,"w_n_",m,".png"))
	fig, ax = subplots(1,1)
	ax.set_xlim([0,2*pi])
	ax.set_ylim([0,2*pi])
	ax.set_xlabel(L"$x_1$", fontsize=30)
	ax.set_ylabel(L"$x_2$", fontsize=30)
	ax.xaxis.set_tick_params(labelsize=30)
	ax.yaxis.set_tick_params(labelsize=30)
	pts = hcat(u...)
	vecs = hcat(v...)
	x_pts = [pts[1,:] - eps*vecs[1,:] pts[1,:] + eps*vecs[1,:]]'
	y_pts = [pts[2,:] - eps*vecs[2,:] pts[2,:] + eps*vecs[2,:]]'

	ax.plot(x_pts, y_pts, "royalblue")
	ax.axis("scaled")
	=#

	return u
end
function tensordot(A, b)
	return [Ai*b[i] for (i, Ai) in enumerate(A)]
end
function pushforward_second_order(u, v1, s)
	d2u = d2step(u, s)
	d2u_v1 = reshape([dot(d2u[:,i], v1) for i=1:4], 2, 2)
	return d2u_v1
end
function create_line_colls(u, q, eps)
	n = length(u)
	lc = zeros(n, 2, 2)
	for i = 1:n
		lc[i,1,:] = u[i] .- eps*q[i]	
		lc[i,2,:] = u[i] .+ eps*q[i]	
	end
	return lc
end
function pushforward(u, q, s)
	du = dstep(u, s)
	q = du*q
	return q
end

