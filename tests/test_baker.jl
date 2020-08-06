include("../examples/baker.jl")
using LinearAlgebra
using JLD
function stable_sens(s4)
d = 2
p = 4

nSteps = 10
u_trj = zeros(d, nSteps)
du_trj = zeros(d,d,nSteps)
s = zeros(p)


vs = zeros(2)
q = rand(2)
q /= norm(q)
nm_q = zeros(nSteps)
le = 0.
# J = cos(4y)

n_exps = size(s4)[1]
dJds = zeros(n_exps)
n_rep = 10
for k=1:n_exps
		s[4] = s4[k]
		
	for j=1:n_rep
		u = 2*pi*rand(d)
		u_trj .= step(u, s, nSteps-1)
		du_trj .= dstep(u_trj, s)
		x = pert(u_trj, 4)

		for i = 1:nSteps
			vs .= du_trj[:,:,i]*vs + x[:,i]
			q .= du_trj[:,:,i]*q 
			
			nm_q[i] = norm(q)
			le += log(nm_q[i])/nSteps
			q ./= nm_q[i]

			vs .-= dot(vs,q)*q
			x1, x2 = u_trj[1,i], u_trj[2,i]
			dJdu = [0., -4*sin(4*x2)]
			dJds[k] += dot(dJdu, vs)/nSteps/n_rep
		end
		#@show le, vs, dJds[k], q
	end
end
return dJds
end
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
		
