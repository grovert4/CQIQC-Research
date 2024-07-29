using YAML, LazyGrids
using MPI
filename = "$(ARGS[1])"
include("$(ARGS[2])")

MPI.Init()
commSize = MPI.Comm_size(MPI.COMM_WORLD)
commRank = MPI.Comm_rank(MPI.COMM_WORLD)
println(filename)
println("Hello from $(commRank) of $(commSize)")
params = YAML.load_file("./Input/$(filename).yml")
#(Uarr, fillingarr) = ndgrid(range(params["U_min"], params["U_max"], params["U_length"]), range(params["filling_min"], params["filling_max"], params["filling_length"]) / (params["filling_max"] * 2))
Us = collect(range(params["U_min"], params["U_max"], params["U_length"]))
#fillings = collect(range(params["filling_min"], params["filling_max"], params["filling_length"]) / (params["filling_max"] * 2))
#Js = collect(range(params["J_min"], params["J_max"], params["J_length"]))
#gridsize = params["J_length"]
# elements_per_process = div(params["J_length"], commSize)
Vs = collect(range(params["V_min"], params["V_max"], params["V_length"]))
gridsize = params["V_length"]
elements_per_process = div(params["V_length"], commSize)
remainder = rem(gridsize, commSize)

start_index = commRank * elements_per_process + min(commRank, remainder) + 1
end_index = start_index + elements_per_process - 1 + (commRank < remainder ? 1 : 0)
Us_r = reverse(Us)
for i in start_index:end_index
    for (j, U) in enumerate(Us_r)
        println(U, " U value")
        params["U"] = round(U, sigdigits=5)
        params["V"] = round(Vs[i], sigdigits=5)
        #params["jh"] = round(Js[i], sigdigits=5)
        #params["filling"] = round(fillings[i], sigdigits=5)
        #params["date"] = filename # change to output file name. 
        if j > 2
            params["U_prev"] = Us_r[j-1]
        else
            delete!(params, "U_prev")
        end
        if params["filling"] == round(0.0, sigdigits=5)
            continue
        else
            MFT(params, filename)
        end
    end
end