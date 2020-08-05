include("../examples/baker.jl")
using PyPlot
using JLD
function obj_fun(x,y)
		return cos(4*y)	
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
function obj_fun_erg_avg(s)
	nSteps = 9999
	u = 2*pi*rand(2)
	J = 0.
	u_trj = step(u,s,nSteps)
	x, y = view(u_trj,1,:), view(u_trj,2,:)
	nSteps = nSteps + 1
	J = sum(obj_fun.(x,y)/nSteps)
	return J
end
function plot_Javg_vs_s(ind)
	fig, ax = subplots(1,1)
	s = zeros(4)
	n_pts = 100
	n_rep = 500
	s2 = LinRange(0.,1.,n_pts)
	J = zeros(n_pts)
	for i = 1:n_pts
		@show s2[i]
		s[ind] = s2[i]
		for n=1:n_rep
			J[i] += obj_fun_erg_avg(s)/n_rep
		end
	end
	save("../data/obj_erg_avg/cos4y_s4.jld",
		 "s$ind", s2,
		"J", J)
	ax.plot(s2, J, ".", ms=4.0)
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax.set_xlabel("\$s_$ind\$",fontsize=28)
	ax.set_ylabel(L"$\langle J\rangle$",fontsize=28)
end

function plot_obj_fun_erg_avg()

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
