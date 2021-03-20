using PyPlot
using JLD
function plot_sens()
	X = load("../data/obj_erg_avg/cos4y_s2.jld")
	s2_arr = X["s2"]
	J_arr = X["J"]
	fig, ax = subplots(1,1)
	ax.plot(s2_arr, J_arr, ".", ms=10.0)
    ax.xaxis.set_tick_params(labelsize=32)
    ax.yaxis.set_tick_params(labelsize=32)
    ax.set_xlabel(L"$s_2$",fontsize=32)
    ax.set_ylabel(L"$\langle J\rangle$",fontsize=32)
	ax.grid(true)
	
	X = load("../data/sens/dJds_s2.jld")

	dJds = X["dJds"]
	s2 = X["s2"]

	X = load("../data/obj_erg_avg/cos4y_s2_sens.jld")
	J = X["J"]

	eps = 2.e-2
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n, 
					2)'
	s_pts = reshape([s2 .- eps s2 .+ eps], n, 
					2)'

	ax.plot(s_pts[:,1:n-1], J_pts[:,1:n-1], "k",lw=6.0)
	ax.plot(s_pts[:,n], J_pts[:,n], "k",lw=6.0,label=L"$d_s\langle J\rangle$ from S3")

	ax.legend(fontsize=32)

	return dJds
	

end
	
