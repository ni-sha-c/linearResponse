using PyPlot
using JLD
function plot_lines()
	X = load("../data/sens/solenoid/dJds_s3_K12.jld")
	rho = X["s"]
	dJds = X["dJds"]

	X = load("../data/obj_erg_avg/solenoid/r2_s3.jld")
    s3 = X["s3"]
    J = X["J"]
    fig, ax = subplots(1,1)
    ax.plot(s3, J, ".", ms=10.0)
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)
    ax.set_xlabel(L"$s_2$",fontsize=28)
    ax.set_ylabel(L"$\langle J\rangle$",fontsize=28)
    ax.grid(true)

    eps = 1.5e-2

	X = load("../data/obj_erg_avg/solenoid/r2_s3_sens.jld")
	J = X["J"]
	s3 = X["s3"]
	
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n,
                    2)'
    s_pts = reshape([s3 .- eps s3 .+ eps], n,
                    2)'

    ax.plot(s_pts[:,1:n-1], J_pts[:,1:n-1], "k",lw=3.0)
    ax.plot(s_pts[:,n], J_pts[:,n], "k",lw=4.0,label=L"$d_s\langle J\rangle$ from S3")

    ax.legend(fontsize=28)


end

