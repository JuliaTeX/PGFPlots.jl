module Plots

export Plot, Histogram, Histogram2, Linear, Linear3, Image, Patch2D, Contour, Scatter, Quiver, Node, Circle, Ellipse, Command

using ..ColorMaps
using Compat
using Discretizers
using StatsBase

typealias RealRange @compat Tuple{Real,Real}

include("ndgrid.jl")

abstract Plot

type Linear <: Plot
    data::AbstractMatrix{Real}
    mark
    markSize
    style
    legendentry
    onlyMarks
    errorBars
    preamble
    epilogue
    Linear{T <: Real}(data::AbstractMatrix{T}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing, errorBars=nothing, preamble=nothing, epilogue=nothing) = new(data, mark, markSize, style, legendentry, onlyMarks, errorBars,preamble,epilogue)
end
Linear{A<:Real, B<:Real}(x::AbstractVector{A}, y::AbstractVector{B}; kwargs...) = Linear(hcat(x, y)'; kwargs...)
Linear{A<:Real}(data::AbstractVector{A}; kwargs...) = Linear(collect(1:length(data)), data; kwargs...)

type Linear3 <: Plot
    data::AbstractMatrix{Real}
    mark
    markSize
    style
    legendentry
    onlyMarks
    preamble
    epilogue
    Linear3{T<:Real}(data::AbstractMatrix{T}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing, preamble=nothing, epilogue=nothing) = new(data, mark, markSize, style, legendentry, onlyMarks,preamble,epilogue)
end
Linear3{A<:Real, B<:Real, C<:Real}(x::AbstractVector{A}, y::AbstractVector{B}, z::AbstractVector{C}; kwargs...) = Linear3(hcat(x, y, z)'; kwargs...)

const THRESHOLD_NSAMPLES_DISC_OURSELVES = 1000 # if we have more samples than this we discretize ourselves
function _construct_histogram_linear_data{Q<:Real,R<:Real}(
    data::Vector{Q},
    binedges::Vector{R},
    density::Bool, # If true, the bar height will be based on the probability density - otherwise directly on counts
    cumulative::Bool, # A cumulative histogram uses the sum of all previous bins and the current one as final value.
    )

    n = length(binedges)
    disc = LinearDiscretizer(binedges)
    counts = get_discretization_counts(disc, data)
    if cumulative
        cumsum!(counts, counts)
    end

    arr_x = convert(Vector{Float64}, binedges)
    arr_y = convert(Vector{Float64}, counts)
    if density
        arr_y ./= sum(counts)
        arr_y ./= binwidths(disc)
    end
    push!(arr_y, arr_y[end])

    Linear(hcat(arr_x, arr_y)', style="ybar interval,fill=blue!10, draw=blue", mark="none")
end
type Histogram <: Plot
    data::AbstractVector{Real}
    bins::Integer
    density::Bool
    cumulative::Bool
    style::AbstractString
    discretization::Symbol
    preamble
    epilogue
    Histogram(data; bins=10, discretization=:default, density=false, cumulative=false, style="fill=blue!10", preamble=nothing, epilogue=nothing) = new(data,bins,density,cumulative,style,discretization,preamble,epilogue)
end

type Contour <: Plot
    data::AbstractMatrix
    xbins
    ybins
    style
    number
    levels
    labels
    preamble
    epilogue
    Contour(data, xbins, ybins; style=nothing, number=nothing, levels=nothing, labels=nothing, preamble=nothing, epilogue=nothing) = new(data, xbins, ybins, style, number, levels, labels,preamble,epilogue)
    function Contour(f::Function, xrange::RealRange, yrange::RealRange; xbins=40, ybins=40, style=nothing, number=nothing, levels=nothing, labels=nothing, preamble=nothing, epilogue=nothing)
        x = linspace(xrange[1], xrange[2], xbins)
        y = linspace(yrange[1], yrange[2], ybins)
        A = zeros(xbins, ybins)
        try
            A = Float64[f(xi, yi) for xi in x, yi in y]
        catch
            A = Float64[f([xi,yi]) for xi in x, yi in y]
        end
        new(A, x, y, style, number, levels, labels, preamble, epilogue)
    end
end

type Scatter <: Plot
    data::AbstractMatrix{Any}
    mark
    markSize
    style
    legendentry
    onlyMarks
    scatterClasses
    preamble
    epilogue
    function Scatter{T<:Any}(data::AbstractMatrix{T}; mark=nothing, markSize=nothing, style=nothing, onlyMarks=true, legendentry=nothing, scatterClasses=nothing, preamble=nothing, epilogue=nothing)
        if size(data,1) == 2
            return Linear(data, mark=mark, markSize=markSize, style=style, onlyMarks=onlyMarks, legendentry=legendentry,preamble=preamble,epilogue=epilogue)
        else
            return new(data, mark, markSize, style, legendentry, onlyMarks, scatterClasses,preamble,epilogue)
        end
    end
end
Scatter{A<:Real, B<:Real}(x::AbstractVector{A}, y::AbstractVector{B}; kwargs...) = Scatter(hcat(x, y)'; kwargs...)
Scatter{A<:Real, B<:Real, C<:Any}(x::AbstractVector{A}, y::AbstractVector{B}, f::AbstractVector{C}; kwargs...) = Scatter(permutedims(hcat(x, y, f), [2,1]); kwargs...)
Scatter{A<:Real, B<:Real}(x::A, y::B; kwargs...) = Scatter(hcat(x, y)'; kwargs...)
Scatter{A<:Real, B<:Real}(x::A, y::B, f; kwargs...) = Scatter(hcat(x, y, f)'; kwargs...)

type Quiver <: Plot
    data::Matrix{Real}
    style
    legendentry
    preamble
    epilogue
    Quiver{T<:Real}(data::Matrix{T}; style=nothing, legendentry=nothing, preamble=nothing, epilogue=nothing) = new(data, style, legendentry,preamble,epilogue)
end
function Quiver(f::Function, xrange::RealRange, yrange::RealRange; style=nothing, legendentry=nothing, samples=15, normalize=true)
    x = linspace(xrange[1], xrange[2], samples)
    y = linspace(yrange[1], yrange[2], samples)
    (X, Y) = meshgrid(x, y)
    n = length(X)
    U = zeros(n)
    V = zeros(n)
    for i = 1:n
        (U[i], V[i]) = f(X[i], Y[i])
    end
    if normalize
        r = max(maximum(U),maximum(V))
        r /= min(minimum(diff(x)),minimum(diff(y)))
        U /= r
        V /= r
    end
    Quiver(X[:], Y[:], U, V, style=style, legendentry=legendentry)
end
Quiver{A<:Real,B<:Real,C<:Real,D<:Real}(x::Vector{A}, y::Vector{B}, u::Vector{C}, v::Vector{D}; kwargs...) = Quiver(hcat(x, y, u, v)'; kwargs...)

type Node <: Plot
    data
    style
    x
    y
    axis # `nothing` will default to "axis cs", other options include "axis description cs", "xticklabel cs", etc.
    preamble
    epilogue
    Node(data, x, y; style=nothing, axis=nothing, preamble=nothing, epilogue=nothing) = new(data, style, x, y, axis,preamble,epilogue)
end

type Circle <: Plot
    xc
    yc
    radius
    style
    preamble
    epilogue
    Circle(xc=0,yc=0,radius=1;style=nothing, preamble=nothing, epilogue=nothing) = new(xc,yc,radius,style,preamble,epilogue)
end

type Ellipse <: Plot
    xc
    yc
    xradius
    yradius
    style
    preamble
    epilogue
    Ellipse(xc=0,yc=0,xradius=1,yradius=1;style=nothing, preamble=nothing, epilogue=nothing) = new(xc,yc,xradius,yradius,style,preamble,epilogue)
end

type Command <: Plot
    cmd::AbstractString
    Command(cmd::AbstractString) = new(cmd)
end

global _imgid = 1

type Image <: Plot
    filename::AbstractString
    xmin::Real
    xmax::Real
    ymin::Real
    ymax::Real
    zmin::Real
    zmax::Real
    colorbar::Bool
    colormap::ColorMaps.ColorMap
    style
    preamble
    epilogue
    function Image{T <: Real}(A::Matrix{T}, xrange::RealRange, yrange::RealRange; filename=nothing, colorbar=true, colormap=ColorMaps.Gray(), zmin=nothing, zmax=nothing, style=nothing, preamble=nothing, epilogue=nothing)
        global _imgid
        if filename == nothing
            id=myid()*10000000000000+_imgid
            filename = "tmp_$(id).png"
            _imgid += 1
        end
        if zmin == nothing
            zmin = minimum(A)
        end
        if zmax == nothing
            zmax = maximum(A)
        end
        if zmin == zmax
            zmin -= 1.
            zmax += 1.
        end
        A = clamp(A, zmin, zmax)
        A = A .- zmin
        A = A ./ (zmax - zmin)
        if isa(colormap, ColorMaps.ColorMap)
            write(colormap, A, filename)
        else
            write(ColorMaps.RGBArray(colormap), A, filename)
        end
        new(filename, xrange[1], xrange[2], yrange[1], yrange[2], zmin, zmax, colorbar, colormap, style, preamble, epilogue)
    end
    function Image(f::Function, xrange::RealRange, yrange::RealRange; filename=nothing, colorbar=true, colormap=ColorMaps.Gray(), zmin=nothing, zmax=nothing, xbins=100, ybins=100, style=nothing, preamble=nothing, epilogue=nothing)
        x = linspace(xrange[1], xrange[2], xbins)
        y = linspace(yrange[1], yrange[2], ybins)
        (X, Y) = meshgrid(x, y)
        A = map(f, X, Y)
        A = flipdim(A, 1)
        Image(A, xrange, yrange, filename=filename, colorbar=colorbar, colormap=colormap, zmin=zmin, zmax=zmax, style=style, preamble=preamble, epilogue=epilogue)
    end
end

type Patch2D{R<:Real} <: Plots.Plot
    data::Matrix{R} # d × m⋅n
                    # where m = 4 for rect and m = 3 for triangle (defaults to triangle)
                    #   and n = number of patches
                    #   and d = {x,y,color} or {x,y}
    style
    patch_type
    shader
    legendentry
end
Patch2D{R<:Real}(data::Matrix{R}; style="patch", patch_type=nothing, shader=nothing, legendentry=nothing) = Patch2D{R}(data, style, patch_type, shader, legendentry)


function Histogram2{A<:Real, B<:Real}(x::Vector{A}, y::Vector{B}; xmin=minimum(x), xmax=maximum(x), ymin=minimum(y), ymax=maximum(y), xbins=50, ybins=50, density=false, filename=nothing, colorbar=true, colormap=ColorMaps.Gray(), zmin=nothing, zmax=nothing, style=nothing)
    ex = linspace(xmin, xmax, xbins+1)
    ey = linspace(ymin, ymax, ybins+1)
    h = fit(StatsBase.Histogram, (y, x), (ey, ex))
    ex, ey, M = h.edges[1], h.edges[2], h.weights
    M = flipdim(M, 1)
    if density
        scale =  xbins * ybins / ((xmax-xmin) * (ymax-ymin) * sum(M))
        M = M * scale
    end
    Image(M, (xmin, xmax), (ymin, ymax), filename=filename, colorbar=colorbar, colormap=colormap, zmin=zmin, zmax=zmax, style=style)
end
function Histogram2{A<:Real, B<:Real, C<:Real}(x::Vector{A}, y::Vector{B}, edges_x::AbstractVector{C}, edges_y::AbstractVector{C};
    density=false,
    style=nothing,
    )

    h = fit(StatsBase.Histogram, (x, y), (edges_x, edges_y))
    ex, ey, M = h.edges[1], h.edges[2], h.weights
    m = length(ex)-1
    n = length(ey)-1
    scale =  m*n / ((ex[end]-ex[1]) * (ey[end]-ey[1]) * sum(M))

    patchdata = Array(Float64, 3, 4*n*m)
    patchidx = 0
    for i in 1 : m
        x₁, x₂ = ex[i], ex[i+1]
        for j in 1 : n
            y₁, y₂ = ey[j], ey[j+1]
            c = M[i,j]
            if density
                c /= (scale*(x₂-x₁)*(y₂-y₁))
            end

            patchidx += 1
            patchdata[1, patchidx] = x₁
            patchdata[2, patchidx] = y₁
            patchdata[3, patchidx] = c

            patchidx += 1
            patchdata[1, patchidx] = x₂
            patchdata[2, patchidx] = y₁
            patchdata[3, patchidx] = c

            patchidx += 1
            patchdata[1, patchidx] = x₂
            patchdata[2, patchidx] = y₂
            patchdata[3, patchidx] = c

            patchidx += 1
            patchdata[1, patchidx] = x₁
            patchdata[2, patchidx] = y₂
            patchdata[3, patchidx] = c
        end
    end

    if isa(style, Void)
        style = "patch"
    elseif isa(style, String)
        style = "patch" * (isempty(style) ? "" : (", "*style))
    end

    Patch2D(patchdata, style=style, patch_type="rectangle")
end

end # end plot module
