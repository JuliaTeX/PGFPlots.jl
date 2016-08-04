module Plots

export Plot, Histogram, Histogram2, Linear, Linear3, ErrorBars, Image, Contour, Scatter, Quiver, Node, Circle, Ellipse, Command

using ..ColorMaps
using Compat
using Discretizers

typealias RealRange @compat Tuple{Real,Real}

include("ndgrid.jl")

abstract Plot

type Linear <: Plot
    data::AbstractArray{Real,2}
    mark
    markSize
    style
    legendentry
    onlyMarks
    Linear{T<:Real}(data::AbstractArray{T,2}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = new(data, mark, markSize, style, legendentry, onlyMarks)
end

type Linear3 <: Plot
    data::AbstractArray{Real,2}
    mark
    markSize
    style
    legendentry
    onlyMarks
    Linear3{T<:Real}(data::AbstractArray{T,2}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = new(data, mark, markSize, style, legendentry, onlyMarks)
end

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
    data::AbstractArray{Real,1}
    bins::Integer
    density::Bool
    cumulative::Bool
    style::AbstractString
end
function Histogram(data; bins=nothing, discretization=:default, density=false, cumulative=false, style="fill=blue!10", mark=nothing, markSize=nothing, legendentry=nothing, onlyMarks=nothing)
    if isa(bins, Integer) && (discretization == :pgfplots == (discretization == :default && length(data) ≤ THRESHOLD_NSAMPLES_DISC_OURSELVES))
        # default - discretize using PGFPlots for smaller sample sizes
        Histogram(data, bins, density, cumulative, style)
    else
        # discretize using Discretizers.jl and produce a Linear plot

        if isa(bins, Integer)
            lo, hi = extrema(data)
            edges = collect(linspace(lo,hi,bins+1))
        else
            if discretization == :default
                discretization = :auto # default is auto
            end
            edges = binedges(DiscretizeUniformWidth(discretization), data)
        end

        linear = _construct_histogram_linear_data(data, edges, density, cumulative)
        if style != "fill=blue!10"
            linear.style = style
            if !contains(linear.style, "ybar interval")
                linear.style = "ybar interval," * linear.style
            end
        end
        if !isa(mark, Void)
            linear.mark = mark
        end
        linear.markSize = markSize
        linear.legendentry = legendentry
        linear.onlyMarks = onlyMarks
        linear
    end
end


type Contour <: Plot
    data::AbstractMatrix{Real}
    xbins
    ybins
    style
    number
    levels
	labels
    Contour(data, xbins, ybins; style=nothing, number=nothing, levels=nothing,labels=nothing) = new(data, xbins, ybins, style, number, levels, labels)
    function Contour(f::Function, xrange::RealRange, yrange::RealRange; xbins=40, ybins=40, style=nothing, number=nothing, levels=nothing, labels=nothing)
        x = linspace(xrange[1], xrange[2], xbins)
        y = linspace(yrange[1], yrange[2], ybins)
        A = zeros(xbins, ybins)
        try
            A = Float64[f(xi, yi) for xi in x, yi in y]
        catch
            A = Float64[f([xi,yi]) for xi in x, yi in y]
        end
        new(A, x, y, style, number, levels, labels)
    end
end

function Contour(z::AbstractMatrix, x::Range, y::Range; style=nothing, number=nothing, levels=nothing, labels=nothing)
    (X, Y) = meshgrid(x, y)
    Contour([X[:]'; Y[:]'; z[:]'], length(x), length(y); style = style, number = number, number = levels, labels=labels)
end

type Scatter <: Plot
    data::AbstractArray{Any,2}
    mark
	markSize
    style
    legendentry
    onlyMarks
	scatterClasses

	function Scatter{T<:Any}(data::AbstractArray{T,2}; mark=nothing, markSize=nothing, style=nothing, onlyMarks=true, legendentry=nothing, scatterClasses=nothing)
		if size(data,1) == 2
			return Linear(data, mark=mark, markSize=markSize, style=style, onlyMarks=onlyMarks, legendentry=legendentry)
		else
			return new(data, mark, markSize, style, legendentry, onlyMarks, scatterClasses)
		end
	end
end

type ErrorBars <: Plot
    data::AbstractArray{Real,2}
    mark
    style
    legendentry
    ErrorBars{T<:Real}(data::AbstractArray{T,2}; mark=nothing, style=nothing, legendentry=nothing) = new(data, mark, style, legendentry)
end

type Quiver <: Plot
    data::Matrix{Real}
    style
    legendentry
    Quiver{T<:Real}(data::Matrix{T}; style=nothing, legendentry=nothing) = new(data, style, legendentry)
end

type Node <: Plot
    data
    style
    x
    y
    axis # `nothing` will default to "axis cs", other options include "axis description cs", "xticklabel cs", etc.
    Node(data, x, y; style=nothing, axis=nothing) = new(data, style, x, y, axis)
end

type Circle <: Plot
	xc
	yc
	radius
    style
	Circle(xc=0,yc=0,radius=1;style=nothing) = new(xc,yc,radius,style)
end

type Ellipse <: Plot
	xc
	yc
	xradius
	yradius
    style
	Ellipse(xc=0,yc=0,xradius=1,yradius=1;style=nothing) = new(xc,yc,xradius,yradius,style)
end

type Command <: Plot
    cmd::AbstractString
    Command(cmd::AbstractString) = new(cmd)
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

Quiver{A<:Real,B<:Real,C<:Real,D<:Real}(x::Vector{A}, y::Vector{B}, u::Vector{C}, v::Vector{D}; style=nothing, legendentry=nothing) = Quiver([x y u v]', style=style, legendentry=legendentry)

Linear{A<:Real, B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = Linear([x y]', mark=mark, markSize=markSize, style=style, legendentry=legendentry, onlyMarks=onlyMarks)
Linear{A<:Real}(data::AbstractArray{A,1}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = Linear([1:length(data)], data, mark=mark, markSize=markSize, style=style, legendentry=legendentry, onlyMarks=onlyMarks)


Linear3{A<:Real, B<:Real, C<:Real}(x::AbstractVector{A}, y::AbstractVector{B}, z::AbstractVector{C}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = Linear3([x y z]', mark=mark, markSize=markSize, style=style, legendentry=legendentry, onlyMarks=onlyMarks)

ErrorBars{A<:Real, B<:Real, C<:Real, D<:Real, E<:Real, F<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}, xplus::AbstractArray{C,1}, yplus::AbstractArray{D,1},
                                                                xminus::AbstractArray{E,1}, yminus::AbstractArray{F,1}; mark=nothing, style=nothing, legendentry=nothing) = ErrorBars([x y xplus yplus xminus yminus]',mark=mark, style=style, legendentry=legendentry)
ErrorBars{A<:Real, B<:Real, C<:Real, D<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}, yplus::AbstractArray{C,1},yminus::AbstractArray{D,1}; mark=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = ErrorBars([x y zeros(length(x)) yplus zeros(length(x)) yminus]', mark=mark, style=style, legendentry=legendentry)
ErrorBars{A<:Real, B<:Real, C<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}, yplusminus::AbstractArray{C,1}; mark=nothing, style=nothing, legendentry=nothing) = ErrorBars([x y zeros(length(x)) yplusminus zeros(length(x)) yplusminus]', mark=mark, style=style, legendentry=legendentry)

Scatter{A<:Real, B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, scatterClasses=nothing) = Scatter([x y]', mark=mark, markSize=markSize, style=style, legendentry=nothing, scatterClasses=scatterClasses)
Scatter{A<:Real, B<:Real, C<:Any}(x::AbstractArray{A,1}, y::AbstractArray{B,1}, f::AbstractArray{C,1}; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, scatterClasses=nothing) = Scatter([x y f]', mark=mark, markSize=markSize, style=style, legendentry=nothing, scatterClasses=scatterClasses)
Scatter{A<:Real, B<:Real}(x::A, y::B; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing) = Scatter([x y]', mark=mark, markSize=markSize, style=style, legendentry=nothing)
Scatter{A<:Real, B<:Real}(x::A, y::B, f; mark=nothing, markSize=nothing, style=nothing, legendentry=nothing, scatterClasses=nothing) = Scatter([x y f]', mark=mark, markSize=markSize, style=style, legendentry=nothing, scatterClasses=scatterClasses)

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
    function Image{T <: Real}(A::Matrix{T}, xrange::RealRange, yrange::RealRange; filename=nothing, colorbar=true, colormap=ColorMaps.Gray(), zmin=nothing, zmax=nothing, style=nothing)
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
        new(filename, xrange[1], xrange[2], yrange[1], yrange[2], zmin, zmax, colorbar, colormap, style)
    end
    function Image(f::Function, xrange::RealRange, yrange::RealRange; filename=nothing, colorbar=true, colormap=ColorMaps.Gray(), zmin=nothing, zmax=nothing, xbins=100, ybins=100, style=nothing)
        x = linspace(xrange[1], xrange[2], xbins)
        y = linspace(yrange[1], yrange[2], ybins)
        (X, Y) = meshgrid(x, y)
        A = map(f, X, Y)
        A = flipdim(A, 1)
        Image(A, xrange, yrange, filename=filename, colorbar=colorbar, colormap=colormap, zmin=zmin, zmax=zmax, style=style)
    end
end

function Histogram2{A<:Real, B<:Real}(x::Vector{A}, y::Vector{B}; xmin=minimum(x), xmax=maximum(x), ymin=minimum(y), ymax=maximum(y), xbins=50, ybins=50, density=false, filename=nothing, colorbar=true, colormap=ColorMaps.Gray(), zmin=nothing, zmax=nothing, style=nothing)
    ex = linspace(xmin, xmax, xbins+1)
    ey = linspace(ymin, ymax, ybins+1)
    ex, ey, M = hist2d(hcat(y, x), ey, ex)
    M = flipdim(M, 1)
    if density
        scale =  xbins * ybins / ((xmax-xmin) * (ymax-ymin) * sum(M))
        M = M * scale
    end
    Image(M, (xmin, xmax), (ymin, ymax), filename=filename, colorbar=colorbar, colormap=colormap, zmin=zmin, zmax=zmax, style=style)
end

end # end plot module
