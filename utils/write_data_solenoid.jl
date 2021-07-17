using JLD
n_min = 32000
n_max = 320000000
filename = "../data/sens/solenoid/dJds_K12_"
function read_data()
    n = n_min
    i = 1
	nn = Integer(log10(n_max/n_min)) + 2
	n_samples = ones(Int64,nn)
	dJds = zeros(nn)
	while n <= n_max
	    X = load(string(filename, n, ".jld"))
	    dJds[i] = X["dJds"][1]
		n_samples[i] = n
	    @show n, dJds[i]
	    i = i + 1
	    n = n*10
    end
	save("../data/sens/solenoid/dJds3_K12_N.jld",
		 "dJds", dJds, "n_samples", n_samples)
end	
