using YAML, LazyGrids
using MPI
include("./Bilayer_MFT_Uniform.jl")
filename = "$(ARGS[1])"
MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)
println(filename)
println("Hello from $(commRank) of $(commSize)")
params = YAML.load_file("./Input/$(filename).yml")

(Uarr, fillingarr) = ndgrid(range(params["U_min"], params["U_max"], params["U_length"]), range(params["filling_min"], params["filling_max"], params["filling_length"]) / 24)
Us = collect(Iterators.flatten(Uarr))
fillings = collect(Iterators.flatten(fillingarr))

gridsize = params["U_length"] * params["filling_length"]

elements_per_process = div(gridsize, commSize)
remainder = rem(gridsize, commSize)

start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)

for i in start_index:end_index
    params["U"] = round(Us[i], sigdigits=5)
    params["filling"] = round(fillings[i], sigdigits=5)
    #params["date"] = filename # change to output file name. 
    MFT(params, filename)
end