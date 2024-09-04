
using TightBindingToolkit, LinearAlgebra
using LaTeXStrings, JLD2
using Plots

using PyPlot
using PyCall

style = pyimport("matplotlib.style")

style.use("lake.mplstyle")
rcParams = pyimport("matplotlib").rcParams

Update the rcParams dictionary
rcParams["text.latex.preamble"] = "\\usepackage{lmodern}"

rcParams["text.usetex"] = true
rcParams["font.family"] = "lmodern"
rcParams["font.serif"] = ["lmodern"] 
rcParams["mathtext.fontset"] = "cm"
rcParams["text.latex.unicode"] = true
PyPlot.rc("text", usetex=true)
function pyplot_hexagonal_old(bz::BZ, data::Matrix{Float64} ;
    colorbar_title::AbstractString=L"\Omega_n(\mathbf{k})",
    cmap::String="Purples", lim::Float64=1.25,
    fontsize::Int=12,
    annotation::AbstractString="",
    annotation_position::Tuple{Float64, Float64}=(0.0, 1.0),
    clims::Tuple{Float64, Float64}=(0.0, 1.0),
    savename::AbstractString="data")

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
    PyPlot.gca().set_xlabel(L"k_x")
    PyPlot.gca().set_ylabel(L"k_y")
    

    if !isempty(savename)
        PyPlot.savefig(savename, bbox_inches="tight")
    end


end
function pyplot_hexagonal(ax::PyCall.PyObject, cax::PyCall.PyObject, input::Dict, data::Matrix{T};
    colorbar_title::AbstractString=L"\Omega_n(\mathbf{k})",
    cmap::String="Blues", lim::Float64=2.025,
    annotation::AbstractString="",
    annotation_position::Tuple{Float64, Float64}=(0.0, 1.0),
    clims::Tuple{Float64, Float64}=(0.0, 5.0),
    savename::AbstractString="") where {T}

    PyPlot.cla()

    im = ax.pcolormesh(input["kxs"], input["kys"], data',
        cmap=cmap, vmin=clims[1], vmax=clims[2])  # Adjust the size as needed

    PyPlot.colorbar(im, label=colorbar_title, cax=cax,fraction=0.1, shrink = 0.5)

    inner_hex = input["inner_hex"]
    for edge in inner_hex
        xs, ys = edge
        ax.plot(xs, ys, linewidth=2.0, c=:orange)
    end

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    PyPlot.annotate(annotation, annotation_position, xycoords="axes fraction")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(L"k_x")
    ax.set_ylabel(L"k_y")

    if !isempty(savename)
        PyPlot.savefig(savename, bbox_inches="tight")
    end


end
params = Dict()
SkXSize = get!(params, "SkXSize", 2)
SkX = get!(params, "SkX", "Bloch")
SkX = "Bloch"
a1 = SkXSize / 2 * [-3.0, sqrt(3)]
a2 = SkXSize / 2 * [3.0, sqrt(3)]
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 2, 2)
for j = 0:(SkXSize-1)
    for i = 0:(SkXSize*3-1)
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end
n = 30
kSize = 6 * n + 3
bz = BZ(kSize)
FillBZ!(bz, UC)

b1 = [bz.basis[1]; 0.0]
b2 = [bz.basis[2]; 0.0]
loc = "/media/andrewhardy/9C33-6BBD/Skyrmion/Monolayer_Data/"
fileName = loc * "J=4.0_band=1_metric.jld2"
fileName_2 = loc * "J=4.0_band=13_metric.jld2"

Data= load(fileName) #MeanFieldToolkit.MFTResume.ReadMFT(fileName)
Data_2= load(fileName_2) #MeanFieldToolkit.MFTResume.ReadMFT(fileName)

kxs = Data["kxs"]
kys = Data["kys"]
curvature = Data["curvature"]
metric_sqrtDets = Data["metric_sqrtDets"]
curvature_2 = Data_2["curvature"]
metric_sqrtDets_2 = Data_2["metric_sqrtDets"]
fig = PyPlot.figure(figsize=(25, 6))
gs = fig.add_gridspec(1, 8, width_ratios=[1,1,1,1,0.15, 0.075, 0.15, 0.075], wspace = 0.01)

ax1 = fig.add_subplot(gs[0, 1])
ax2 = fig.add_subplot(gs[0, 2])
ax3 = fig.add_subplot(gs[0, 3])
ax4 = fig.add_subplot(gs[0, 4])
cax1 = fig.add_subplot(gs[0, 6])
cax2 = fig.add_subplot(gs[0, 8])

pyplot_hexagonal(ax1, cax1, Data, abs.(curvature), cmap = "Purples")
pyplot_hexagonal(ax2,cax2, Data, metric_sqrtDets, cmap = "Blues")
pyplot_hexagonal(ax3,cax1, Data_2, abs.(curvature_2), cmap = "Purples")
pyplot_hexagonal(ax4,cax2, Data_2, metric_sqrtDets_2, cmap = "Blues")

display(gcf())

