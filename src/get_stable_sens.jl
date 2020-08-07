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
	
