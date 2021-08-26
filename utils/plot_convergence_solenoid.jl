using PyPlot
using JLD
nSteps = 200001
X = load(string("../data/sens/solenoid/metrics_",nSteps,"_K.jld"))
dJds_ust = X["dJds_ust"]
dJds_st = X["dJds_st"]
ust_cn_k = X["unstable_contribution_k"]
K = size(ust_cn_k)[1]
fig, ax = subplots()
ust_pge = abs.((cumsum(ust_cn_k) .- dJds_ust)./dJds_ust)*100 
ax.semilogy(1:K, ust_pge, "ko-", ms=20, lw=2.5, label="% error in sum of k terms")
ax.set_xlabel("k",fontsize=28)
ax.xaxis.set_tick_params(labelsize=28)
ax.yaxis.set_tick_params(labelsize=28)
ax.grid(true)
leg = ax.legend(fontsize=28)
leg.legendHandles[1]._legmarker.set_markersize(10)



