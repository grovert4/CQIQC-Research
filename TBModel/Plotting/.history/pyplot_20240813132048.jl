function pyplot_hexagonal_old(bz::BZ, data::Matrix{Float64} ;
    colorbar_title::AbstractString=L"\Omega_n(\mathbf{k})",
    cmap::String="Purples", lim::Float64=1.25,
    fontsize::Int=12,
    annotation::AbstractString="",
    annotation_position::Tuple{Float64, Float64}=(0.0, 1.0),
    clims::Tuple{Float64, Float64}=(0.0, 1.0),
    savename::AbstractString="")

    PyPlot.cla()
    normalized_data = data ./ 1.0
    b1, b2 = bz.basis
    kSize = bz.gridSize[1]

    shifts = (kSize-1)/kSize .* [[0.0, 0.0], b1, -b1, b2, -b2, b1+b2, -b1-b2, b1-b2, b2-b1]

    for shift in shifts
        PyPlot.pcolormesh(getindex.(bz.ks .+ Ref(shift), 1), getindex.(bz.ks .+ Ref(shift), 2), normalized_data,
            cmap=cmap, vmin=clims[1], vmax=clims[2])  # Adjust the size as needed
    end
    PyPlot.colorbar(label=colorbar_title)

    RotMat = [cos(pi/3) sin(pi/3); -sin(pi/3) cos(pi/3)]
    corner = bz.HighSymPoints["K1"]

    for i in 1:6
        p1 = corner
        p2 = RotMat * corner

        if abs(p2[1]-p1[1])>1e-6
            slope = (p2[2] - p1[2]) / (p2[1] - p1[1])
            xs = LinRange(p1[1], p2[1], 100)
            ys = slope .* (xs .- p1[1]) .+ p1[2]
            PyPlot.plot(xs, ys, linewidth=2.0, c=:black)
        else
            ys = LinRange(p1[2], p2[2], 100)
            xs = p1[1] .* ones(100)
            PyPlot.plot(xs, ys, linewidth=2.0, c=:black)
        end
        corner = RotMat * corner
    end

    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().set_xticks(collect(LinRange(-5, 5, 11)))
    PyPlot.gca().set_yticks(collect(LinRange(-5, 5, 11)))
    PyPlot.annotate(annotation, annotation_position, xycoords="axes fraction")
    PyPlot.gca().set_xlim(-lim, lim)
    PyPlot.gca().set_ylim(-lim, lim)
    PyPlot.gca().set_xlabel(L"k_x", fontfamily = "Helvetica")
    PyPlot.gca().set_ylabel(L"k_y", fontfamily = "Helvetica")

    if !isempty(savename)
        PyPlot.savefig(savename, bbox_inches="tight")
    end


end
function pyplot_hexagonal(input::Dict, data::Matrix{T};
    colorbar_title::AbstractString=L"\Omega_n(\mathbf{k})",
    cmap::String="Blues", lim::Float64=2.025,
    annotation::AbstractString="",
    annotation_position::Tuple{Float64, Float64}=(0.0, 1.0),
    clims::Tuple{Float64, Float64}=(0.0, 7.0),
    savename::AbstractString="") where {T}

    PyPlot.cla()

    PyPlot.pcolormesh(input["kxs"], input["kys"], data',
        cmap=cmap, vmin=clims[1], vmax=clims[2])  # Adjust the size as needed

    PyPlot.colorbar(label=colorbar_title)

    inner_hex = input["inner_hex"]
    for edge in inner_hex
        xs, ys = edge
        PyPlot.plot(xs, ys, linewidth=2.0, c=:orange)
    end

    PyPlot.gca().set_aspect("equal")
    PyPlot.gca().set_xticks([])
    PyPlot.gca().set_yticks([])
    PyPlot.annotate(annotation, annotation_position, xycoords="axes fraction")
    PyPlot.gca().set_xlim(-lim, lim)
    PyPlot.gca().set_ylim(-lim, lim)
    PyPlot.gca().set_xlabel(L"k_x")
    PyPlot.gca().set_ylabel(L"k_y")

    if !isempty(savename)
        PyPlot.savefig(savename, bbox_inches="tight")
    end


end