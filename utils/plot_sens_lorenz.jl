using PyPlot
using JLD
function plot_errorbar()
	X = load("../data/sens/lorenz/dJds_vs_params.jld")
	rho_arr = X["s"]
	dJds_arr = X["dJds"]
	err_arr = sqrt.(X["var_dJds"])/4
	@show err_arr
	fig, ax = subplots(1,1)
	ax.plot(rho_arr, dJds_arr, "X", ms=10.0)
	ax.errorbar(rho_arr, dJds_arr, yerr=err_arr, linewidth=3.0, color="blue",fmt="none")
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)
    ax.set_xlabel(L"$s$",fontsize=28)
	ax.set_ylabel(L"$d_{s}\langle J\rangle$",fontsize=28)
	ax.grid(true)
	

end
function plot_lines()
	X = load("../data/sens/lorenz/dvarz_ds2_800.jld")
	rho = X["rho"]
	dJds = X["dJds"]
	err = sqrt.(X["var_dJds"])/4

	dJds = dJds[rho <= 31]
	err = err[rho <= 31]


	@show err
	X = load("../data/obj_erg_avg/lorenz/zm28sq.jld")
    s2 = X["s2"]
    J = X["J"]
    fig, ax = subplots(1,1)
    ax.plot(s2, J, ".", ms=20.0)
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)
    ax.set_xlabel(L"$s$",fontsize=28)
    ax.set_ylabel(L"$\langle J\rangle$",fontsize=28)
    ax.grid(true)

    eps = 1.2e-1

	X = load("../data/obj_erg_avg/lorenz/zm28sq_dJds.jld")
	J = X["J"]
	s2 = X["s2"]

	J = J[rho <= 31]
	s2 = s2[rho <= 31]
	rho = rho[rho <= 31]
	s2 = s2[err .<= 1.8]
	J = J[err .<= 1.8]
    dJds = dJds[err .<= 1.8]
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n,
                    2)'
    s_pts = reshape([s2 .- eps s2 .+ eps], n,
                    2)'

    ax.plot(s_pts[:,1:n-1], J_pts[:,1:n-1], "k",lw=3.0)
    ax.plot(s_pts[:,n], J_pts[:,n], "k",lw=3.0,label=L"$d_s\langle J\rangle$ from S3")

    ax.legend(fontsize=28)

	X = load("../data/sens/lorenz/dvarz_ds2_800.jld")
	rho = X["rho"]
	dJds = X["dJds"]
	err = sqrt.(X["var_dJds"])/4
	@show err
	X = load("../data/obj_erg_avg/lorenz/zm28sq_dJds.jld")
	J = X["J"]
	s2 = X["s2"]

	s2 = s2[err .<= 1.8]
	J = J[err .<= 1.8]
    dJds = dJds[err .<= 1.8]
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n,
                    2)'
    s_pts = reshape([s2 .- eps s2 .+ eps], n,
                    2)'

    ax.plot(s_pts[:,1:n-1], J_pts[:,1:n-1], "k",lw=3.0)
    ax.plot(s_pts[:,n], J_pts[:,n], "k",lw=3.0,label=L"$d_s\langle J\rangle$ from S3")




end

