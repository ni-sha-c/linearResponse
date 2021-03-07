using PyPlot
using JLD
function plot_stable_sens()
	X = load("../data/obj_erg_avg/cos4y_s4.jld")
	s4_arr = X["s4"]
	J_arr = X["J"]
	fig, ax = subplots(1,1)
	ax.plot(s4_arr, J_arr, ".", ms=10.0)
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)
    ax.set_xlabel(L"$s_4$",fontsize=28)
    ax.set_ylabel(L"$\langle J\rangle$",fontsize=28)
	ax.grid(true)
	
	X = load("../data/stable_sens/dJds4.jld")
	dJds = X["dJds"]
	s4 = X["s4"]

	X = load("../data/obj_erg_avg/cos4y_s4_sens.jld")
	J = X["J"]

	eps = 2.e-2
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n, 
					2)'
	s_pts = reshape([s4 .- eps s4 .+ eps], n, 
					2)'

	ax.plot(s_pts[:,1:3:end-1], J_pts[:,1:3:end-1], "k",lw=6.0)
	ax.plot(s_pts[:,end],J_pts[:,end],"k",lw=6.0,label=L"$d_s\langle J\rangle$ from S3")
	ax.legend(fontsize=28)

end
	
