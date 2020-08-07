using JLD
function plot_stable_sens()
	X = load("../data/obj_erg_avg/sin2x_cos4y_s4.jld")
	s4_arr = X["s4"]
	J_arr = X["J"]
	fig, ax = subplots(1,1)
	ax.plot(s4_arr, J_arr, ".", ms=10.0)
    ax.xaxis.set_tick_params(labelsize=28)
    ax.yaxis.set_tick_params(labelsize=28)
    ax.set_xlabel(L"$s_4$",fontsize=28)
    ax.set_ylabel(L"$\langle J\rangle$",fontsize=28)
	ax.grid(true)
	
	n_fac = 10
	s4 = s4_arr[1:n_fac:end]
	dJds = stable_sens(s4)
	J = J_arr[1:n_fac:end]
	eps = 1.e-2
	n = size(dJds)[1]
	J_pts = reshape([J .- eps*dJds J .+ eps*dJds], n, 
					2)'
	s_pts = reshape([s4 .- eps s4 .+ eps], n, 
					2)'

	ax.plot(s_pts, J_pts, "g",lw=2.0)



end
	
