using PyPlot
using JLD
function plot_sens()
	X = load("../data/obj_erg_avg/cos4y_s1_s3.jld")
	s1_arr = X["s1"]
	J_arr = X["J"]
	fig, ax = subplots(1,1)
	ax.plot(s1_arr, J_arr, ".", ms=10.0)
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)
    ax.set_xlabel(L"$s_1, s_3$",fontsize=28)
    ax.set_ylabel(L"$\langle J\rangle$",fontsize=28)
	ax.grid(true)
	
	X = load("../data/sens/dJds.jld")

	dJds = X["dJds"]
	s1 = X["s1"]

	X = load("../data/obj_erg_avg/cos4y_s1_s3_sens.jld")
	J = X["J"]

	eps = 1.e-2
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n, 
					2)'
	s_pts = reshape([s1 .- eps s1 .+ eps], n, 
					2)'

	ax.plot(s_pts[:,1:n-1], J_pts[:,1:n-1], "k",lw=6.0)
	ax.plot(s_pts[:,n], J_pts[:,n], "k",lw=6.0,label=L"$d_s\langle J\rangle$ from S3")

	ax.legend(fontsize=28)

	return dJds
	

end
	
