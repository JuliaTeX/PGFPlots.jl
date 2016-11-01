VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module PGFPlots

export LaTeXString, @L_str, @L_mstr
export plot, ErrorBars, Axis, Axes, PolarAxis, GroupPlot, Plots, ColorMaps, save, define_color
export pushPGFPlotsOptions, popPGFPlotsOptions, resetPGFPlotsOptions, pgfplotsoptions
export pushPGFPlotsPreamble, popPGFPlotsPreamble, resetPGFPlotsPreamble, pgfplotspreamble
export pushPGFPlots, popPGFPlots
import Colors: RGB
import Contour: contours, levels
import Reel: extension

using Compat
using Discretizers
using DataFrames

type ErrorBars
    data::AbstractMatrix{Real}
    style
    mark
    ErrorBars{T <: Real}(data::AbstractMatrix{T}; style=nothing, mark=nothing) = new(data, style, mark)
end

mylength(x) = (x == nothing) ? 0 : length(x)

function ErrorBars(; x=nothing, y=nothing, xplus=nothing, xminus=nothing, yplus=nothing, yminus=nothing, style=nothing, mark=nothing)
    if x != nothing
        xplus = x
        xminus = x
    end
    if y != nothing
        yplus = y
        yminus = y
    end
    n = maximum(map(mylength, Any[xplus, xminus, yplus, yminus]))
    xplus = mylength(xplus) == n ? xplus : zeros(n)
    xminus = mylength(xminus) == n ? xminus : zeros(n)
    yplus = mylength(yplus) == n ? yplus : zeros(n)
    yminus = mylength(yminus) == n ? yminus : zeros(n)
    data = hcat(xplus, xminus, yplus, yminus)'
    ErrorBars(data, style=style, mark=mark)
end

include("colormaps.jl")
include("plots.jl")

Plots.Linear(df::DataFrame; x::Symbol=:x, y::Symbol=:y, kwargs...) = Plots.Linear(df[x], df[y]; kwargs...)
Plots.Linear3(df::DataFrame; x::Symbol=:x, y::Symbol=:y, z::Symbol=:z, kwargs...) = Plots.Linear3(df[x], df[y], df[z]; kwargs...)
Plots.Histogram(df::DataFrame; x::Symbol=:x, kwargs...) = Plots.Histogram(df[:x]; kwargs...)
Plots.Scatter(df::DataFrame; x::Symbol=:x, y::Symbol=:y, kwargs...) = Plots.Scatter(df[x], df[y]; kwargs...)
Plots.Quiver(df::DataFrame; x::Symbol=:x, y::Symbol=:y, u::Symbol=:u, v::Symbol=:v, kwargs...) = Plots.Quiver(convert(Vector, df[x]), convert(Vector, df[y]), convert(Vector, df[u]), convert(Vector, df[v]); kwargs...)
Plots.Histogram2(df::DataFrame; x::Symbol=:x, y::Symbol=:y, kwargs...) = Plots.Histogram2(convert(Vector, df[x]), convert(Vector, df[y]); kwargs...)
Plots.Histogram2{C<:Real}(df::DataFrame, edges_x::AbstractVector{C}, edges_y::AbstractVector{C}; x::Symbol=:x, y::Symbol=:y, kwargs...) = Plots.Histogram2(convert(Vector, df[x]), convert(Vector, df[y]), edges_x, edges_y; kwargs...)

import TikzPictures: TikzPicture, PDF, TEX, SVG, save, LaTeXString, @L_str, @L_mstr

_pgfplotsoptions = Any[""]

pushPGFPlotsOptions(options::AbstractString) = push!(_pgfplotsoptions, options)
function popPGFPlotsOptions()
    if length(_pgfplotsoptions) > 1
        pop!(_pgfplotsoptions)
    end
end
resetPGFPlotsOptions() = deleteat!(_pgfplotsoptions, 2:length(_pgfplotsoptions))

function pushPGFPlots(options::AbstractString)
    warn("Use pushPGFPlotsOptions")
    pushPGFPlotsOptions(options)
end

function popPGFPlots()
    warn("Use popPGFPlotsOptions")
    popPGFPlotsOptions()
end

pgfplotsoptions() = join(_pgfplotsoptions, ",\n")

_pgfplotspreamble = Any[readstring(joinpath(dirname(@__FILE__), "preamble.tex"))]

pushPGFPlotsPreamble(preamble::AbstractString) = push!(_pgfplotspreamble, preamble)
function popPGFPlotsPreamble()
    if length(_pgfplotspreamble) > 1
        pop!(_pgfplotspreamble)
    end
end
resetPGFPlotsPreamble() = deleteat!(_pgfplotspreamble,2:length(_pgfplotspreamble))

pgfplotspreamble() = join(_pgfplotspreamble, "\n")

function define_color(name::AbstractString, color::RGB)
    _pgfplotspreamble[end] = _pgfplotspreamble[end] * "\n\\definecolor{$name}{rgb}{$(color.r), $(color.g), $(color.b)}"
end

define_color(name::AbstractString, color) = define_color(name, convert(RGB, color))

function define_color{T<:AbstractFloat}(name::AbstractString, color::Vector{T})
    if length(color) == 3
        r, g, b = color[1], color[2], color[3]
        _pgfplotspreamble[end] = _pgfplotspreamble[end] * "\n\\definecolor{$name}{rgb}{$r, $g, $b}"
    elseif length(color) == 4
        c, m, y, k = color[1], color[2], color[3], color[4]
        _pgfplotspreamble[end] = _pgfplotspreamble[end] * "\n\\definecolor{$name}{cmyk}{$c, $m, $y, $k}"
    else
        error("Color must have three or four components")
    end
end

function define_color{T<:Integer}(name::AbstractString, color::Vector{T})
    @assert(length(color) == 3, "Color must have three components")
    r, g, b = color[1], color[2], color[3]
    _pgfplotspreamble[end] = _pgfplotspreamble[end] * "\n\\definecolor{$name}{RGB}{$r, $g, $b}"
end

function define_color(name::AbstractString, color::UInt32)
    _pgfplotspreamble[end] = _pgfplotspreamble[end] * "\n\\definecolor{$name}{HTML}{$(uppercase(hex(color)))}"
end

function define_color(name::AbstractString, color::Real)
    _pgfplotspreamble[end] = _pgfplotspreamble[end] * "\n\\definecolor{$name}{gray}{$color}"
end

histogramMap = Dict(
    :bins => "bins",
    :density => "density",
    :cumulative => "cumulative"
    )

linearMap = Dict(
    :mark => "mark",
    :markSize => "mark size",
    :style => "",
    :onlyMarks => "only marks"
    )

scatterMap = Dict(
    :mark => "mark",
    :markSize => "mark size",
    :style => "",
    :onlyMarks => "only marks",
    :scatterClasses => "scatter/classes"
    )

errorbarsMap = Dict(
    :mark => "error mark",
    :style => "error bar style"
    )

quiverMap = Dict(
    :style => ""
    )

contourMap = Dict(
    :number => "number",
    :levels => "levels",
    :style => "",
    :labels => "labels"
    )

patch2DMap = Dict(
    :style => "",
    :patch_type => "patch type",
    :shader => "shader",
    :legendentry => "legendentry"
    )

using .Plots
using .ColorMaps

typealias IntegerRange @compat Tuple{Integer,Integer}
typealias RealRange @compat Tuple{Real,Real}

type Axis
    plots::Vector{Plot}
    title
    xlabel
    xlabelStyle
    ylabel
    ylabelStyle
    zlabel
    zlabelStyle
    xmin
    xmax
    ymin
    ymax
    axisEqual
    axisEqualImage
    enlargelimits
    axisOnTop
    view
    width
    height
    style
    legendPos
    legendStyle
    xmode
    ymode
    colorbar
    hideAxis
    axisLines
    axisKeyword

    Axis{P <: Plot}(plots::Vector{P};title=nothing, xlabel=nothing, xlabelStyle=nothing, ylabel=nothing, ylabelStyle=nothing, zlabel=nothing, zlabelStyle=nothing, xmin=nothing, xmax=nothing,
                    ymin=nothing, ymax=nothing, axisEqual=nothing, axisEqualImage=nothing, enlargelimits=nothing, axisOnTop=nothing, view=nothing, width=nothing,
                    height=nothing, style=nothing, legendPos=nothing, legendStyle=nothing, xmode=nothing, ymode=nothing, colorbar=nothing, hideAxis=nothing, axisLines=nothing, axisKeyword="axis") =
        new(plots, title, xlabel, xlabelStyle, ylabel, ylabelStyle, zlabel, zlabelStyle, xmin, xmax, ymin, ymax, axisEqual, axisEqualImage, enlargelimits, axisOnTop, view, width, height, style, legendPos, legendStyle, xmode, ymode, colorbar, hideAxis, axisLines, axisKeyword
            )

    Axis(;kwargs...) = Axis(Plot[]; kwargs...)
    Axis(plot::Plot; kwargs...) = Axis(Plot[plot]; kwargs...)

    Axis(plot::Contour; kwargs...) = Axis(Plot[plot]; kwargs..., view="{0}{90}")
    Axis(plots::Vector{Contour}; kwargs...) = Axis(Plot[plots...]; kwargs..., view="{0}{90}")
end

PolarAxis(args...; kwargs...) = Axis(args...; kwargs..., axisKeyword = "polaraxis")

typealias Axes Vector{Axis}

function Base.push!(g::Axis, p::Plot)
    push!(g.plots, p)
end

axisMap = Dict(
    :title => "title",
    :xlabel => "xlabel",
    :xlabelStyle => "xlabel style",
    :ylabel => "ylabel",
    :ylabelStyle => "ylabel style",
    :zlabel => "zlabel",
    :zlabelStyle => "zlabel style",
    :xmin => "xmin",
    :xmax => "xmax",
    :ymin => "ymin",
    :ymax => "ymax",
    :axisEqual => "axis equal",
    :axisEqualImage => "axis equal image",
    :enlargelimits => "enlargelimits",
    :axisOnTop => "axis on top",
    :view => "view",
    :width => "width",
    :height => "height",
    :style => "",
    :legendPos => "legend pos",
    :legendStyle => "legend style",
    :xmode => "xmode",
    :ymode => "ymode",
    :colorbar => "colorbar",
    :hideAxis => "hide axis",
    :axisLines => "axis lines"
    )

type GroupPlot
    axes::AbstractArray{Axis,1}
    dimensions::IntegerRange
    style
    groupStyle
    GroupPlot(;style=nothing,groupStyle=nothing) = new(AbstractArray{Axis,1}[], (0,0), style, groupStyle)
    GroupPlot(axes::AbstractArray{Axis,1};style=nothing,groupStyle=nothing) = new(axes, (length(axes), 1), style, groupStyle)
    GroupPlot(plots;style=nothing,groupStyle=nothing) = new([Axis(p) for p in plots], (length(plots), 1), style, groupStyle)
    GroupPlot(dimensions::IntegerRange;style=nothing,groupStyle=nothing) = new(AbstractArray{Axis,1}[], dimensions, style, groupStyle)
    GroupPlot(rows::Integer, columns::Integer;style=nothing,groupStyle=nothing) = new(AbstractArray{Axis,1}[], (rows, columns), style, groupStyle)
end

function Base.push!(g::GroupPlot, a::Axis)
    push!(g.axes, a)
    if length(g.axes) > prod(g.dimensions)
        g.dimensions = (length(g.axes), 1)
    end
    g
end
Base.push!(g::GroupPlot, p::Plot) = push!(g, Axis(p))
function Base.append!{T<:Union{Plot,Axis}}(g::GroupPlot, arr::AbstractVector{T})
    for a in arr
        push!(g, a)
    end
    g
end

function printList{T}(o::IOBuffer, a::AbstractArray{T,1}; brackets=false)
    first = true
    for elem in a
        if first
            first = false
            if brackets
                print(o, "{")
            end
            print(o, "$elem")
        else
            print(o, ", $elem")
        end
    end
    if !first && brackets
        print(o, "}")
    end
end

printObject(o::IOBuffer, object) = print(o, "$(object)")
printObject{T}(o::IOBuffer, object::AbstractArray{T,1}) = printList(o, object, brackets = true)

function optionHelper(o::IOBuffer, m, object; brackets=false, otherOptions=Dict{AbstractString,AbstractString}[], otherText=nothing)
    first = true
    for (sym, str) in m
        if getfield(object,sym) != nothing
            if first
                first = false
                if brackets
                    print(o, "[")
                end
            else
                print(o, ", ")
            end
            if length(str) > 0
                print(o, "$str = {")
            end
            printObject(o, getfield(object,sym))
            if length(str) > 0
                print(o, "}")
            end
        end
    end
    for (k, v) in otherOptions
        if first
            first = false
            if brackets
                print(o, "[")
            end
        else
            print(o, ", ")
        end
        print(o, "$k = $v")
    end
    if otherText != nothing
        for t in otherText
            if t != nothing
                if first
                    first = false
                    if brackets
                        print(o, "[")
                    end
                else
                    print(o, ", ")
                end
                print(o, "$t")
            end
        end
    end
    if !first && brackets
        print(o, "]")
    end
end

function plotHelper(o::IOBuffer, p::Plots.Histogram)
    if (p.discretization == :pgfplots && p.bins > 0) ||
       (p.discretization == :default && length(p.data) ≤ Plots.THRESHOLD_NSAMPLES_DISC_OURSELVES)

        print(o, "\\addplot+ [mark=none, $(p.style), hist={")
        optionHelper(o, histogramMap, p)
        print(o, "}] table [row sep=\\\\, y index = 0] {")
        println(o, "data\\\\")
        for d in p.data
            println(o, "$d \\\\ ")
        end
        println(o, "};")

    else
        # discretize using Discretizers.jl

        if p.discretization == :specified && p.bins > 0
            edges = binedges(DiscretizeUniformWidth(p.bins), convert(Vector{Float64}, p.data))
        else
            discretization = p.discretization == :default ? :auto : p.discretization
            edges = binedges(DiscretizeUniformWidth(discretization), convert(Vector{Float64}, p.data))
        end

        linear = Plots._construct_histogram_linear_data(p.data, edges, p.density, p.cumulative)
        if p.style != "fill=blue!10"
            linear.style = p.style
            if !contains(linear.style, "ybar interval")
                linear.style = "ybar interval," * linear.style
            end
        end
        plotHelper(o, linear)
    end
end

plotLegend(o::IOBuffer, entry) = nothing
plotLegend(o::IOBuffer, entry::AbstractString) = println(o, "\\addlegendentry{$entry}")
function plotLegend{T <: AbstractString}(o::IOBuffer, entries::Vector{T})
    for entry in entries
        plotLegend(o, entry)
    end
end

# todo: add error bars style
function plotHelperErrorBars(o::IOBuffer, p::Linear)
    println(o, "[")
    optionHelper(o, linearMap, p, brackets=false)
    println(o, ", error bars/.cd, ")
    optionHelper(o, errorbarsMap, p.errorBars, otherText=["x dir=both", "x explicit", "y dir=both", "y explicit"])
    println(o, "]")
    x = vcat(p.data, p.errorBars.data)
    println(o, "table [")
    println(o, "x error plus=ex+, x error minus=ex-, y error plus=ey+, y error minus=ey-")
    println(o, "] {")
    println(o, "x y ex+ ex- ey+ ey-")
    for i = 1:size(p.data, 2)
        println(o, "$(x[1,i]) $(x[2,i]) $(x[3,i]) $(x[4,i]) $(x[5,i]) $(x[6,i])")
    end
    println(o, "};")
end

function plotHelper(o::IOBuffer, p::Linear)
    print(o, p.preamble)
    print(o, "\\addplot+ ")
    if p.errorBars == nothing
        optionHelper(o, linearMap, p, brackets=true)
        println(o, "coordinates {")
        for i = 1:size(p.data,2)
            println(o, "($(p.data[1,i]), $(p.data[2,i]))")
        end
        println(o, "};")
    else
        plotHelperErrorBars(o, p)
    end
    plotLegend(o, p.legendentry)
    print(o, p.epilogue)
end

function plotHelper(o::IOBuffer, p::Scatter)
    print(o, p.preamble)
    if p.scatterClasses == nothing
        print(o, "\\addplot+[scatter, scatter src=explicit, ")
    else
        print(o, "\\addplot+[scatter, scatter src=explicit symbolic, ")
    end
    optionHelper(o, scatterMap, p)
    println(o, "] coordinates {")
    for i = 1:size(p.data,2)
        println(o, "($(p.data[1,i]), $(p.data[2,i])) [$(p.data[3,i])]")
    end
    println(o, "};")
    plotLegend(o, p.legendentry)
    print(o, p.epilogue)
end

# Specific version for Linear3 type
# Changes are addplot3 (vs addplot) and iterate over all 3 columns
function plotHelper(o::IOBuffer, p::Linear3)
    print(o, p.preamble)
    print(o, "\\addplot3+ ")
    optionHelper(o, linearMap, p, brackets=true)
    println(o, "coordinates {")
    for i = 1:size(p.data,2)
        println(o, "($(p.data[1,i]), $(p.data[2,i]), $(p.data[3,i]))")
    end
    println(o, "};")
    plotLegend(o, p.legendentry)
    print(o, p.epilogue)
end

function plotHelper(o::IOBuffer, p::Node)
    print(o, p.preamble)
    axis = p.axis != nothing ? p.axis : "axis cs"

    if p.style != nothing
        println(o, "\\node at ($(axis):$(p.x), $(p.y)) [$(p.style)] {$(p.data)};")
    else
        println(o, "\\node at ($(axis):$(p.x), $(p.y)) {$(p.data)};")
    end
    print(o, p.epilogue)
end


function plotHelper(o::IOBuffer, p::ErrorBars)
    print(o, p.preamble)
    print(o, "\\addplot+ [")
    optionHelper(o, errorbarsMap, p)
    print(o, ",error bars/.cd, x dir=both, x explicit, y dir=both, y explicit")
    println(o,"]")
    println(o, "coordinates {")
    for i = 1:size(p.data,2)
        println(o, "($(p.data[1,i]), $(p.data[2,i])) +=($(p.data[3,i]),$(p.data[4,i])) -=($(p.data[5,i]),$(p.data[6,i]))")
    end
    println(o, "};")
    plotLegend(o, p.legendentry)
    print(o, p.epilogue)
end

function plotHelper(o::IOBuffer, p::Quiver)
    print(o, p.preamble)
    print(o, "\\addplot+ ")
    optionHelper(o, quiverMap, p, brackets=true, otherOptions=Dict("quiver"=>"{u=\\thisrow{u},v=\\thisrow{v}}"))
    println(o, "table {")
    println(o, "x y u v")
    for i = 1:size(p.data,2)
        println(o, "$(p.data[1,i]) $(p.data[2,i]) $(p.data[3,i]) $(p.data[4,i])")
    end
    println(o, "};")
    plotLegend(o, p.legendentry)
    print(o, p.epilogue)
end

function plotHelper(o::IOBuffer, p::Contour)
    print(o, p.preamble)
    arg = 5
    if p.number != nothing
        arg = p.number
    elseif p.levels != nothing
        arg = p.levels
    end
    C = contours(convert(Vector{Float64}, p.xbins), convert(Vector{Float64}, p.ybins), convert(Matrix{Float64}, p.data), arg)
    if p.labels == nothing
        p.labels = true
    end
    if p.labels
        if p.style != nothing
            print(o, "\\addplot3[contour prepared, $(p.style)] table {")
        else
            print(o, "\\addplot3[contour prepared] table {")
        end
    else
        if p.style != nothing
            print(o, "\\addplot3[contour prepared={labels=false}, $(p.style)] table {")
        else
            print(o, "\\addplot3[contour prepared={labels=false}] table {")
        end
    end
    for c in levels(C)
        level = c.level
        for l in c.lines
            for v in l.vertices
                println(o, "$(v[1]) $(v[2]) $level")
            end
            println(o)
        end
    end
    println(o, "};")
    print(o, p.epilogue)
end

function plotHelper(o::IOBuffer, p::Circle)
    print(o, p.preamble)
    if p.style != nothing
        println(o, "\\draw[$(p.style)] (axis cs:$(p.xc), $(p.yc)) circle[radius=$(p.radius)];")
    else
        println(o, "\\draw (axis cs:$(p.xc), $(p.yc)) circle[radius=$(p.radius)];")
    end
    print(o, p.epilogue)
end

function plotHelper(o::IOBuffer, p::Ellipse)
    print(o, p.preamble)
    if p.style != nothing
        println(o, "\\draw[$(p.style)] (axis cs:$(p.xc), $(p.yc)) ellipse[x radius=$(p.xradius), y radius=$(p.yradius)];")
    else
        println(o, "\\draw (axis cs:$(p.xc), $(p.yc)) ellipse[x radius=$(p.xradius), y radius=$(p.yradius)];")
    end
    print(o, p.epilogue)
end

function plotHelper(o::IOBuffer, p::Command)
    println(o, p.cmd*";") #PGFPlots expects commands to be terminated with a ;
end

function plotHelper(o::IOBuffer, p::Image)
    print(o, p.preamble)
    if p.zmin == p.zmax
        error("Your colorbar range limits must not be equal to each other.")
    end
    if p.style != nothing
        println(o, "\\addplot [$(p.style), point meta min=$(p.zmin), point meta max=$(p.zmax)] graphics [xmin=$(p.xmin), xmax=$(p.xmax), ymin=$(p.ymin), ymax=$(p.ymax)] {$(p.filename)};")
    else
        println(o, "\\addplot [point meta min=$(p.zmin), point meta max=$(p.zmax)] graphics [xmin=$(p.xmin), xmax=$(p.xmax), ymin=$(p.ymin), ymax=$(p.ymax)] {$(p.filename)};")
    end
    print(o, p.epilogue)
end

function plotHelper{R<:Real}(o::IOBuffer, p::Patch2D{R})
    print(o, "\\addplot")

    optionHelper(o, patch2DMap, p, brackets=true)
    println(o)

    m = p.patch_type == "rectangle" ? 4 : 3

    if R <: AbstractFloat
        println(o, "table[point meta=\\thisrow{c}] {")
        println(o, "\tx y c")
        for i in 1 : size(p.data, 2)
            @printf(o, "\t %.6e %.6e", convert(Float64, p.data[1,i]), # x
                                       convert(Float64, p.data[2,i])) # y
            if size(p.data, 1) > 2
                @printf(o, " %.6e\n", convert(Float64, p.data[3,i])) # color
            else
                println(o, "")
            end

            if mod(i, m) == 0
                println(o) # extra space between patches (trainges or rectangles)
            end
        end
    else
        println(o, "table {")
        println(o, "\tx y")
        for i in 1 : size(p.data, 2)
            @printf(o, "\t %d %d", p.data[1,i], # x
                                   p.data[2,i]) # y
            if size(p.data, 1) > 2
                @printf(o, " %d\n", p.data[3,i]) # color
            else
                println(o, "")
            end

            if mod(i, m) == 0
                println(o) # extra space between patches (trainges or rectangles)
            end
        end
    end
    println(o, "};")
    plotLegend(o, p.legendentry)
end

# plot option string and contents; no \begin{axis} or \nextgroupplot
function plotHelper(o::IOBuffer, axis::Axis)
    optionHelper(o, axisMap, axis, brackets=true, otherText=[axisOptions(p) for p in axis.plots])
    for p in axis.plots
        plotHelper(o, p)
    end
end

function plot(axis::Axis)
    o = IOBuffer()
    print(o, "\\begin{$(axis.axisKeyword)}")
    plotHelper(o, axis)
    println(o, "\\end{$(axis.axisKeyword)}")
    TikzPicture(takebuf_string(o), options=pgfplotsoptions(), preamble=pgfplotspreamble())
end

function plot(axes::Axes)
    o = IOBuffer()

    for axis in axes
        print(o, "\\begin{$(axis.axisKeyword)}")
        plotHelper(o, axis)
        println(o, "\\end{$(axis.axisKeyword)}")
    end
    TikzPicture(takebuf_string(o), options=pgfplotsoptions(), preamble=pgfplotspreamble())
end

function plot(p::GroupPlot)
    o = IOBuffer()
    style = ""
    if p.style != nothing
        style = p.style * ", "
    end
    groupStyle = ""
    if p.groupStyle != nothing
        groupStyle = p.groupStyle * ", "
    end
    println(o, "\\begin{groupplot}[$(style)group style={$(groupStyle)group size=$(p.dimensions[1]) by $(p.dimensions[2])}]")
    for a in p.axes
        print(o, "\\nextgroupplot ")
        plotHelper(o, a)
    end
    println(o, "\\end{groupplot}")
    mypreamble = pgfplotspreamble() * "\\usepgfplotslibrary{groupplots}"
    TikzPicture(takebuf_string(o), options=pgfplotsoptions(), preamble=mypreamble)
end

typealias Plottable Union{Plot,GroupPlot,Axis,Axes,TikzPicture}

plot(p::Plot) = plot(Axis(p))

function plot(p::Contour)
    plot(Axis(p, xmin=p.xbins[1], xmax=p.xbins[end], ymin=p.ybins[1], ymax=p.ybins[end]))
end

plot{A<:Real,B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}; kwargs...) = plot(Linear(x, y; kwargs...))

plot{A<:Real,B<:Real,C<:Real}(x::AbstractVector{A}, y::AbstractVector{B}, z::AbstractVector{C}; kwargs...) = plot(Linear3(x, y, z; kwargs...))

function Plots.Linear(f::Function, range::RealRange; mark="none", style=nothing, legendentry=nothing)
    x = linspace(range[1], range[2])
    y = map(f, x)
    Linear(x, y, mark=mark, style=style, legendentry=legendentry)
end

plot(f::Function, range::RealRange; kwargs...) = plot(Linear(f, range; kwargs...))

plot(tkz::TikzPicture) = tkz # tikz pic doesn't need plot, here for convenience

Base.mimewritable(::MIME"image/svg+xml", p::Plottable) = true

cleanup(p::Axis) = map(cleanup, p.plots)
cleanup(axes::Axes) = map(cleanup, axes)
cleanup(p::GroupPlot) = map(cleanup, p.axes)
cleanup(p::Plot) = nothing
cleanup(p::Circle) = nothing
cleanup(p::Ellipse) = nothing
cleanup(p::Command) = nothing
cleanup(p::Image) = rm(p.filename)
cleanup(p::Contour) = nothing

axisOptions(p::Plot) = nothing

function colormapOptions(cm::ColorMaps.Gray)
    if cm.invert
        return "colormap={wb}{gray(0cm)=(1); gray(1cm)=(0)}"
    else
        return "colormap={wb}{gray(0cm)=(0); gray(1cm)=(1)}"
    end
end

function colormapOptions(cm::ColorMaps.RGBArray)
    o = IOBuffer()
    print(o, "colormap={mycolormap}{ ")
    n = length(cm.colors)
    for i = 1:n
        c = cm.colors[i]
        print(o, "rgb($(i-1)cm)=($(c.r),$(c.g),$(c.b)) ")
    end
    print(o, "}")
    takebuf_string(o)
end

function axisOptions(p::Image)
    if p.colorbar
        cmOpt = colormapOptions(p.colormap)
        return "enlargelimits = false, axis on top, $cmOpt, colorbar"
    else
        return "enlargelimits = false, axis on top"
    end
end

@compat function Base.show(f::IO, a::MIME"image/svg+xml", p::Plottable)
    r = Base.show(f, a, plot(p))
    cleanup(p)
    r
end

function save(filename::AbstractString, o::Plottable; include_preamble::Bool=true)
    _, ext = splitext(filename)
    ext = lowercase(ext)
    if ext == ".pdf"
        save(PDF(filename), plot(o))
    elseif ext == ".svg"
        save(SVG(filename), plot(o))
    elseif ext == ".tex"
        save(TEX(filename, include_preamble=include_preamble), plot(o))
    elseif ext == "." || ext == ""
        error("You must specify a file extension.")
    else
        error("Unsupported file extensions: $ext")
    end
end

###################################################
# Reel configuration

extension(m::MIME"image/svg+xml") = "svg"

end # module
