using PyPlot
using JLD
nSteps = 100001
X = load("../data/sens/solenoid/metrics_$nSteps.jld")
dJds_ust = X["dJds_ust"]
dJds_st = X["dJds_st"]
x_trj = X["x_trj"]
t_trj = (atan.(x_trj[2,:],x_trj[1,:]) .+ 
		 2π) .% 2π
g = X["g"]
a = X["a"]
b = X["b"]
nv = X["normvs"]

runup = 5000
t_trj = t_trj[runup:end]
g = g[runup:end]
a = a[runup:end]
b = b[runup:end]
nv = nv[runup:end]


fig, ax = subplots()
ax.plot(t_trj, g, "r.", ms=0.5, label=L"g")
ax.axis("scaled")
ax.xaxis.set_tick_params(labelsize=28)
ax.yaxis.set_tick_params(labelsize=28)
ax.grid(true)
leg = ax.legend(fontsize=28)
leg.legendHandles[1]._legmarker.set_markersize(10)


fig, ax = subplots()
ax.plot(t_trj, a, "r.", ms=0.5, label=L"a")
ax.axis("scaled")
ax.xaxis.set_tick_params(labelsize=28)
ax.yaxis.set_tick_params(labelsize=28)
ax.grid(true)
leg = ax.legend(fontsize=28)
leg.legendHandles[1]._legmarker.set_markersize(10)


fig, ax = subplots()
ax.plot(t_trj, b, "r.", ms=0.5, label=L"b")
ax.axis("scaled")
ax.xaxis.set_tick_params(labelsize=28)
ax.yaxis.set_tick_params(labelsize=28)
ax.grid(true)
ax.set_xlabel(L"\theta", fontsize=28)
leg = ax.legend(fontsize=28)
leg.legendHandles[1]._legmarker.set_markersize(10)



fig, ax = subplots()
ax.plot(t_trj, nv, "b.", ms=0.5, label=L"\||v\||")
ax.set_xlabel(L"\theta", fontsize=28)
ax.xaxis.set_tick_params(labelsize=28)
ax.yaxis.set_tick_params(labelsize=28)
ax.grid(true)
leg = ax.legend(fontsize=28)
leg.legendHandles[1]._legmarker.set_markersize(10)


