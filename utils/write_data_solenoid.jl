using JLD
dirname = "../data/sens/solenoid/dJds_vs_N"
function read_data()
    filenames = readdir(dirname)
	npts = size(filenames)[1]
	n_samples = zeros(Int64, npts)
	dJds = zeros(npts)
	for (i, filename) in enumerate(filenames)
			n_samples[i] = parse(Int64, split(filename, ".")[1])
			X = load(string(dirname,"/",filename))
			dJds[i] = X["dJds"][1]
	end

	save("../data/sens/solenoid/dJds3.jld",
		 "dJds", dJds, "n_samples", n_samples)
end	
