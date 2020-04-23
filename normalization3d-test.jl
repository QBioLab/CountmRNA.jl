using FileIO
include("normalization3d.jl")
include("tiffxml.jl")


data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200316_result"
#cell_name = filter(x->occursin(r"^s",x), readdir(data_dir) )
cell_name = filter(x->occursin(r"^s[0-9][0-9]?-c[0-9][0-9]?.tiff-clear.ome.tiff$",x), readdir(data_dir) )

parameters = []
for name in cell_name
	print("processing $name")
	@time cell = load(File(format"TIFF", "$data_dir/$name"))
	@time imgs_norms, para = normalize(cell)
	push!(parameters, para)
end

save("normalization3d-test.jld", "parameters", parameters)
