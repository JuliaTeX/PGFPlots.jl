module Plots

export Plot, Histogram, Histogram2, BarChart, Linear, Linear3, Image, Patch2D, Contour, Scatter, Quiver, Node, Circle, Ellipse, Command

using ..ColorMaps
using Discretizers
using StatsBase
using Distributed

const RealRange = Tuple{Real,Real}

include("ndgrid.jl")

abstract type Plot end

mutable struct Linear <: Plot
    data::AbstractMatrix{Real}
    mark
    markSize
    style
    legendentry
    onlyMarks
    errorBars
    Linear(data::AbstractMatrix{T}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing, errorBars=nothing) where {T <: Real} = new(data, mark, markSize, style, legendentry, onlyMarks, errorBars)
end
Linear(x::AbstractVector{A}, y::AbstractVector{B}; kwargs...) where {A<:Real, B<:Real} = Linear(hcat(x, y)'; kwargs...)
Linear(data::AbstractVector{A}; kwargs...) where {A<:Real} = Linear(collect(1:length(data)), data; kwargs...)

mutable struct Linear3 <: Plot
    data::AbstractMatrix{Real}
    mark
    markSize
    style
    legendentry
    onlyMarks
    Linear3(data::AbstractMatrix{T}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) where {T<:Real} = new(data, mark, markSize, style, legendentry, onlyMarks)
end
Linear3(x::AbstractVector{A}, y::AbstractVector{B}, z::AbstractVector{C}; kwargs...) where {A<:Real, B<:Real, C<:Real} = Linear3(hcat(x, y, z)'; kwargs...)

const THRESHOLD_NSAMPLES_DISC_OURSELVES = 1000 # if we have more samples than this we discretize ourselves
function _construct_histogram_linear_data(
    data::Vector{Q},
    binedges::Vector{R},
    density::Bool, # If true, the bar height will be based on the probability density - otherwise directly on counts
    cumulative::Bool, # A cumulative histogram uses the sum of all previous bins and the current one as final value.
    ) where {Q<:Real,R<:Real}

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
mutable struct Histogram <: Plot
    data::AbstractVector{Real}
    bins::Integer
    density::Bool
    cumulative::Bool
    style::AbstractString
    discretization::Symbol
    Histogram(data; bins=10, discretization=:default, density=false, cumulative=false, style="fill=blue!10") = new(data,bins,density,cumulative,style,discretization)
end

mutable struct BarChart <: Plot
    keys # symbolic x coords
    values
    style
    legendentry
    errorBars

    BarChart(keys::AbstractVector, values::AbstractVector{R}; style=nothing, legendentry=nothing, errorBars=nothing) where {R<:Real} = new(keys, values, style, legendentry, errorBars)
end
function BarChart(
    values::AbstractVector{R};
    kwargs...) where {R<:Real}

    keys = [string(i) for i in 1 : length(values)]
    return BarChart(keys, values; kwargs...)
end
function BarChart(
    values::AbstractVector{S};
    kwargs...) where {S<:AbstractString}

    dict = Dict{S,Int}()
    for v in values
        dict[v] = get(dict, v, 0) + 1
    end

    keys = S[]
    values = Int[]
    for (k,v) in dict
        push!(keys, k)
        push!(values, v)
    end

    return BarChart(keys, values; kwargs...)
end
function BarChart(
    values::AbstractVector{V},
    disc::CategoricalDiscretizer{N,D}; kwargs...) where {N,D,V}

    subkeys = collect(keys(disc.n2d))
    subvalues = encode(disc, values)
    return BarChart(keys, values; kwargs...)
end
function symbolic_x_coords(keys::AbstractVector)
    retval = ""
    for (i,k) in enumerate(keys)
        retval *= string(k)
        if i != length(keys)
            retval *= ", "
        end
    end
    return retval
end
symbolic_x_coords(p::BarChart) = symbolic_x_coords(p.keys)

mutable struct Contour <: Plot
    data::AbstractMatrix
    xbins
    ybins
    style
    contour_style
    number
    levels
    labels
    Contour(data, xbins, ybins; style=nothing, contour_style=nothing, number=nothing, levels=nothing, labels=nothing) = new(data, xbins, ybins, style, contour_style, number, levels, labels)
    function Contour(f::Function, xrange::RealRange, yrange::RealRange; xbins=40, ybins=40, style=nothing, contour_style=nothing, number=nothing, levels=nothing, labels=nothing)
        x = range(xrange[1], stop=xrange[2], length=xbins)
        y = range(yrange[1], stop=yrange[2], length=ybins)
        A = zeros(xbins, ybins)
        try
            A = Float64[f(xi, yi) for xi in x, yi in y]
        catch
            A = Float64[f([xi,yi]) for xi in x, yi in y]
        end
        new(A, x, y, style, contour_style, number, levels, labels)
    end
end

mutable struct Scatter <: Plot
    data::AbstractMatrix{Any}
    mark
    markSize
    style
    legendentry
    onlyMarks
    scatterClasses

    function Scatter(data::AbstractMatrix; mark=nothing, markSize=nothing, style=nothing, onlyMarks=true, legendentry=nothing, scatterClasses=nothing)
        new(data, mark, markSize, style, legendentry, onlyMarks, scatterClasses)
    end
end
Scatter(x::AbstractVector{A}, y::AbstractVector{B}; kwargs...) where {A<:Real, B<:Real} = Scatter(hcat(x, y)'; kwargs...)
Scatter(x::AbstractVector{A}, y::AbstractVector{B}, f::AbstractVector{C}; kwargs...) where {A<:Real, B<:Real, C<:Any} = Scatter(permutedims(hcat(x, y, f), [2,1]); kwargs...)
Scatter(x::A, y::B; kwargs...) where {A<:Real, B<:Real} = Scatter(hcat(x, y)'; kwargs...)
Scatter(x::A, y::B, f; kwargs...) where {A<:Real, B<:Real} = Scatter(hcat(x, y, f)'; kwargs...)

mutable struct Quiver <: Plot
    data::AbstractMatrix{Real}
    style
    legendentry
    Quiver(data::AbstractMatrix{T}; style=nothing, legendentry=nothing) where {T<:Real} = new(data, style, legendentry)
end
function Quiver(f::Function, xrange::RealRange, yrange::RealRange; style=nothing, legendentry=nothing, samples=15, normalize=true)
    x = range(xrange[1], stop=xrange[2], length=samples)
    y = range(yrange[1], stop=yrange[2], length=samples)
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
function Quiver(
    x::Vector{A},
    y::Vector{B},
    u::Vector{C},
    v::Vector{D};
    kwargs...
    ) where {A<:Real,B<:Real,C<:Real,D<:Real}
    return Quiver(hcat(x, y, u, v)'; kwargs...)
end

mutable struct Node <: Plot
    data
    style
    x
    y
    axis # `nothing` will default to "axis cs", other options include "axis description cs", "xticklabel cs", etc.
    Node(data, x, y; style=nothing, axis=nothing) = new(data, style, x, y, axis)
end

mutable struct Circle <: Plot
    xc
    yc
    radius
    style
    Circle(xc=0,yc=0,radius=1;style=nothing) = new(xc,yc,radius,style)
end

mutable struct Ellipse <: Plot
    xc
    yc
    xradius
    yradius
    style
    Ellipse(xc=0,yc=0,xradius=1,yradius=1;style=nothing) = new(xc,yc,xradius,yradius,style)
end

mutable struct Command <: Plot
    cmd::AbstractString
    Command(cmd::AbstractString) = new(cmd)
end

global _imgid = 1

mutable struct Image <: Plot
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
    function Image(
        A::AbstractMatrix{T},
        xrange::RealRange,
        yrange::RealRange;
        filename=nothing,
        colorbar=true,
        colormap=ColorMaps.GrayMap(),
        zmin=nothing,
        zmax=nothing,
        style=nothing,
        ) where {T <: Real}

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
        A = clamp.(A, zmin, zmax)
        A = A .- zmin
        A = A ./ (zmax - zmin)
        if isa(colormap, ColorMaps.ColorMap)
            write(colormap, A, filename)
        else
            write(ColorMaps.RGBArrayMap(colormap), A, filename)
        end
        new(filename, xrange[1], xrange[2], yrange[1], yrange[2], zmin, zmax, colorbar, colormap, style)
    end
    function Image(f::Function, xrange::RealRange, yrange::RealRange; filename=nothing, colorbar=true, colormap=ColorMaps.GrayMap(), zmin=nothing, zmax=nothing, xbins=100, ybins=100, style=nothing)
        x = range(xrange[1], stop=xrange[2], length=xbins)
        y = range(yrange[1], stop=yrange[2], length=ybins)
        (X, Y) = meshgrid(x, y)
        A = map(f, X, Y)
        A = reverse(A, dims=1)
        Image(A, xrange, yrange, filename=filename, colorbar=colorbar, colormap=colormap, zmin=zmin, zmax=zmax, style=style)
    end
end

mutable struct Patch2D <: Plots.Plot
    data::AbstractMatrix{Real} # d × m⋅n
                    # where m = 4 for rect and m = 3 for triangle (defaults to triangle)
                    #   and n = number of patches
                    #   and d = {x,y,color} or {x,y}
    style
    patch_type
    shader
    legendentry
end
Patch2D(data::AbstractMatrix; style="patch", patch_type=nothing, shader=nothing, legendentry=nothing) = Patch2D(data, style, patch_type, shader, legendentry)


function Histogram2(
    x::Vector{A},
    y::Vector{B};
    xmin=minimum(x),
    xmax=maximum(x),
    ymin=minimum(y),
    ymax=maximum(y),
    xbins=50,
    ybins=50,
    density=false,
    zmode=nothing,
    filename=nothing,
    colorbar=true,
    colormap=ColorMaps.GrayMap(),
    zmin=nothing,
    zmax=nothing,
    style=nothing,
    ) where {A<:Real, B<:Real}

    ex = range(xmin, stop=xmax, length=xbins+1)
    ey = range(ymin, stop=ymax, length=ybins+1)
    h = fit(StatsBase.Histogram, (y, x), (ey, ex), closed=:left)
    ex, ey, M = h.edges[1], h.edges[2], h.weights
    M = reverse(M, dims=1)
    if density
        scale =  xbins * ybins / ((xmax-xmin) * (ymax-ymin) * sum(M))
        M = M * scale
    end
    if zmode == "log" 
        M = M .+ 1
        M = log10.(M)
    end
    Image(M, (xmin, xmax), (ymin, ymax), filename=filename, colorbar=colorbar, colormap=colormap, zmin=zmin, zmax=zmax, style=style)
end
function Histogram2(
    x::Vector{A},
    y::Vector{B},
    edges_x::AbstractVector{C},
    edges_y::AbstractVector{C};
    density=false,
    style=nothing,
    ) where {A<:Real, B<:Real, C<:Real}

    h = fit(StatsBase.Histogram, (x, y), (edges_x, edges_y), closed=:left)
    ex, ey, M = h.edges[1], h.edges[2], h.weights
    m = length(ex)-1
    n = length(ey)-1
    scale =  m*n / ((ex[end]-ex[1]) * (ey[end]-ey[1]) * sum(M))

    patchdata = Array{Float64}(undef, 3, 4*n*m)
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

    if isa(style, Nothing)
        style = "patch"
    elseif isa(style, String)
        style = "patch" * (isempty(style) ? "" : (", "*style))
    end

    Patch2D(patchdata, style=style, patch_type="rectangle")
end

end # end plot module
