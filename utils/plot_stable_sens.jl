using PyPlot
using JLD
function plot_stable_sens()
	X = load("../data/obj_erg_avg/cos4y_s1.jld")
	s4_arr = X["s1"]
	J_arr = X["J"]
	fig, ax = subplots(1,1)
	ax.plot(s4_arr, J_arr, "x", ms=4.0)
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)
    ax.set_xlabel(L"$s_1$",fontsize=28)
    ax.set_ylabel(L"$\langle J\rangle$",fontsize=28)
	ax.grid(true)
	
	X = load("../data/unstable_sens/dJds1.jld")
	dJds = X["dJds"]
	s4 = X["s1"]

	X = load("../data/obj_erg_avg/cos4y_s1_sens.jld")
	J = X["J"]

	eps = 3.e-2
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n, 
					2)'
	s_pts = reshape([s4 .- eps s4 .+ eps], n, 
					2)'

	ax.plot(s_pts[:,1:5:end], J_pts[:,1:5:end], "g",lw=2.0)



end
	