using PyPlot
using JLD
function plot_sens_centerstable()
	X = load("../data/sens/lorenz/dJds_unstable.jld")
	rho_arr = X["rho"]
	dJds_arr = X["dJds"]
	err_arr = sqrt.(X["var_dJds"])/4
	@show err_arr
	fig, ax = subplots(1,1)
	ax.plot(rho_arr, dJds_arr, ".", ms=10.0)
    ax.xaxis.set_tick_params(labelsize=32)
    ax.yaxis.set_tick_params(labelsize=32)
    ax.set_xlabel(L"$s_2$",fontsize=32)
	ax.set_ylabel(L"$d_{s_2}\langle J\rangle$",fontsize=32)
	ax.grid(true)
	

end
	
