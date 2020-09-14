include("../examples/baker.jl")
using LinearAlgebra
using PyPlot
using PyCall
@pyimport matplotlib.collections as coll 
@pyimport matplotlib.colors as cs 
@pyimport matplotlib.cm as cm 
function plot_y(m)
	s = [0.1, 0.1, 0., 0.1]
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
	eps1 = 1.e-3
	filepath = string("/home/nishac/Research/PhDThesis/papers/",
					  "Decomposition-of-Linear-Response/figs/")
	for i = 1:20
		@show i
		u .= next.(u, Ref(s))
		q .= pushforward.(u, q, Ref(s))
		z .= norm.(q)
		q .= q./z
	end
	y_pos = 6.0
	x = Array{Float64}[]
	yx = Array{Array{Float64}}[]
	vx = Array{Array{Float64}}[]


	for i = 1:m
		@show i
		D2 .= pushforward_second_order.(u, q, Ref(s))	
		y .= pushforward.(u, y, Ref(s))
		y .+= tensordot(D2, v)
		X .= pert.(u, 1) .+ pert.(u, 2) .+ pert.(u, 4) 
		DX .= dpert.(u,1) .+ dpert.(u, 2) .+ dpert.(u, 4)
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
		if i >=5
			inds = [abs(uk[2] - y_pos) < eps1 for uk in u]
			push!(x, [ui[1] for ui in u[inds]])
			push!(yx, y[inds])
			push!(vx, v[inds])
		end

	end
	qperp = [[-qi[2], qi[1]] for qi in q]
	ny = dot.(y, qperp) 
	x = vcat(x...)
	yx = vcat(yx...)
	vx = vcat(vx...)
	vxy = [vxi[2] for vxi in vx] 
	yxy = [yxi[2] for yxi in yx] 
	fig, ax = subplots(1,1)
	ax.set_xlim([0,2*pi])
	ax.set_xlabel(L"$x_1$", fontsize=30)
	ax.xaxis.set_tick_params(labelsize=30)
	ax.yaxis.set_tick_params(labelsize=30)
	plt.tight_layout()
	ax.plot(x, vxy,".", label=L"$v_{x_1}(x_2 = 6.0)$", ms=2.0, color="royalblue")
	#ax.plot([(x .- eps)'[1:10:end] .% (2*pi); 
	#		(x .+ eps)'[1:10:end] .% (2*pi)], 
	#		[(vxy .- eps*yxy)'[1:10:end]; 
	#		 (vxy .+ eps*yxy)'[1:10:end]],
	#		color="blue")  
	#ax.plot(x, vx,".", color="darkcyan")

	ax.set_aspect(1.5)
	ax.grid(true)
	leg = ax.legend(loc="right", bbox_to_anchor=(0.8,0.2), 
					fontsize=30, framealpha=0)
	leg.legendHandles[1]._legmarker.set_markersize(20)
	return x, vxy, yx
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
	eps1 = 1.e-3
	filepath = string("/home/nishac/Research/PhDThesis/papers/",
					  "Decomposition-of-Linear-Response/figs/")
	y_pos = 4.5
	x = Array{Float64}[]
	qx = Array{Array{Float64}}[]
	vx = Array{Array{Float64}}[]
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
		if i >=5
			inds = [abs(uk[2] - y_pos) < eps1 for uk in u]
			push!(x, [ui[1] for ui in u[inds]])
			push!(qx, q[inds])
			push!(vx, v[inds])
		end

	end
	x = vcat(x...)
	qx = vcat(qx...)
	vx = vcat(vx...)
		
	fig, ax = subplots(1,1)
	ax.set_xlim([0,2*pi])
	ax.set_xlabel(L"$x_1$", fontsize=30)
	ax.xaxis.set_tick_params(labelsize=30)
	ax.yaxis.set_tick_params(labelsize=30)
	plt.tight_layout()
	ax.plot(x, vx,".",ms=0.5,label=(L"$v_{x_1}$", L"$v_{x_2}$"))
	#ax.plot(x, vx,".", color="darkcyan")
	fig.legend(fontsize=30)
	#ax.axis("scaled")

	return x, vx
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

