using PyPlot
function plot_smooth_indicator()
	fig = figure() 
	ax = fig.add_subplot(111, projection="3d")
	eps = 1.e-8
	n_xpts, n_ypts = 100,100
	x1 = LinRange(0., 2*pi-eps, n_xpts)
	x2 = LinRange(0., 2*pi-eps, n_ypts)
	h1, h2 = 1.0, 1.0
	x1_gr = repeat(x1, 1, n_ypts)
	x2_gr = repeat(x2', n_xpts, 1)
	n_ctrs = 10
	a = 2*pi*rand(n_ctrs)
	b = 2*pi*rand(n_ctrs)
	x3_gr = zeros(n_xpts, n_ypts, n_ctrs)
	for (i, a1) = enumerate(a)
		x3_gr[:,:,i] = smooth_indicator.(x1_gr,
						x2_gr, a1, b[i], h1, h2)
		ax.plot_surface(x1_gr, x2_gr, x3_gr[:,:,i], 
				label="\$ a = $a1 \$",cmap="Blues",
				vmin=0.,vmax=1.0)
	end
	ax.xaxis.set_tick_params(labelsize=28)
	ax.yaxis.set_tick_params(labelsize=28)
	ax.zaxis.set_tick_params(labelsize=28)
	ax.set_xlabel("x",fontsize=28)
	ax.set_ylabel("y",fontsize=28)
	ax.set_zlabel("z",fontsize=28)
	ax.set_title("\$ a = $(a[1]) \$",fontsize=28)
	#lgnd = fig.legend(fontsize=28)
end
function smooth_indicator(x1,x2,a1,a2,h1,h2)
	f1, f2 = 0., 0.
	if (-h1/2 <= x1 - a1 <= h1/2)
		f1 = cos(pi*(x1-a1)/h1)
	end
	if	(-h1/2 <= (x1 - 2*pi) - a1 <= h1/2) 
		f1 = cos(pi*(x1-a1-2*pi)/h1)
	end
	if (-h1/2 <= (x1 + 2*pi) - a1 <= h1/2) 
		f1 = cos(pi*(x1-a1+2*pi)/h1)
	end
	if (-h2/2 <= x2 - a2 <= h2/2) 
		f2 = cos(pi*(x2-a2)/h2)
	end
	if	(-h2/2 <= (x2 - 2*pi) - a2 <= h2/2) 
		f2 = cos(pi*(x2-a2-2*pi)/h2)
	end
	if	(-h2/2 <= (x2 + 2*pi) - a2 <= h2/2)
		f2 = cos(pi*(x2-a2+2*pi)/h2)
	end
	return f1*f2
end
