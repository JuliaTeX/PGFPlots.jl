module PGFPlots

export LaTeXString, @L_str, @L_mstr
export plot, Axis, PolarAxis, GroupPlot, Plots, save, pushPGFPlots, popPGFPlots

include("plots.jl")

import TikzPictures: TikzPicture, PDF, TEX, SVG, save, LaTeXString, @L_str, @L_mstr

_pgfplotsoptions = [""]

pushPGFPlots(options::String) = push!(_pgfplotsoptions, options)
function popPGFPlots()
  if length(_pgfplotsoptions) > 1
    pop!(_pgfplotsoptions)
  end
end

pgfplotsoptions() = _pgfplotsoptions[end]

function pgfplotsoptions(options::String)
  global _pgfplotsoptions
  _pgfplotsoptions = options
end

preamble = readall(joinpath(Pkg.dir("PGFPlots"), "src", "preamble.tex"))

histogramMap = [
  :bins => "bins",
  :density => "density",
  :cumulative => "cumulative"
  ]

linearMap = [
  :mark => "mark",
  :style => "",
  :onlyMarks => "only marks"
  ]

errorbarsMap = [
  :mark => "mark",
  :style => ""
  ]

quiverMap = [
  :style => ""
  ]


contourMap = [
  :number => "number",
  :levels => "levels",
  :style => ""
]


using .Plots

type Axis
  plots::Vector{Plot}
  title
  xlabel
  ylabel
  xmin
  xmax
  ymin
  ymax
  enlargelimits
  axisOnTop
  view
  width
  height
  style
  legendPos

  Axis(plot::Plot;title=nothing, xlabel=nothing, ylabel=nothing, xmin=nothing, xmax=nothing,
       ymin=nothing, ymax=nothing, enlargelimits=nothing, axisOnTop=nothing, view="{0}{90}", width=nothing, height=nothing, style=nothing, legendPos=nothing) =
    new([plot], title, xlabel, ylabel, xmin, xmax, ymin, ymax, enlargelimits, axisOnTop, view, width, height, style, legendPos
  )
  Axis{P <: Plot}(plots::Vector{P};title=nothing, xlabel=nothing, ylabel=nothing, xmin=nothing, xmax=nothing,
       ymin=nothing, ymax=nothing, enlargelimits=nothing, axisOnTop=nothing, view="{0}{90}", width=nothing, height=nothing, style=nothing, legendPos=nothing) =
    new(plots, title, xlabel, ylabel, xmin, xmax, ymin, ymax, enlargelimits, axisOnTop, view, width, height, style, legendPos
  )

end

type PolarAxis
  plots::Vector{Plot}
  title
  yticklabel

  PolarAxis(plot::Plot;title=nothing, yticklabel=nothing) = new([plot], title, yticklabel)
  PolarAxis{P <: Plot}(plots::Vector{P};title=nothing, yticklabel=nothing) = new(plots, title, yticklabel)
end

axisMap = [
  :title => "title",
  :xlabel => "xlabel",
  :ylabel => "ylabel",
  :xmin => "xmin",
  :xmax => "xmax",
  :ymin => "ymin",
  :ymax => "ymax",
  :enlargelimits => "enlargelimits",
  :axisOnTop => "axis on top",
  :view => "view",
  :width => "width",
  :height => "height",
  :style => "",
  :legendPos => "legend pos"
  ]

polarAxisMap = [
  :title => "title",
  :yticklabel => "yticklabel"
  ]

type GroupPlot
  axes::AbstractArray{Axis,1}
  dimensions::(Integer,Integer)
  GroupPlot() = new(AbstractArray{Axis,1}[], (0,0))
  GroupPlot(axes::AbstractArray{Axis,1}) = new(axes, (length(axes), 1))
  GroupPlot(plots) = new([Axis(p) for p in plots], (length(plots), 1))
  GroupPlot(dimensions::(Integer,Integer)) = new(AbstractArray{Axis,1}[], dimensions)
  GroupPlot(rows::Integer, columns::Integer) = new(AbstractArray{Axis,1}[], (rows, columns))
end

function Base.push!(g::GroupPlot, a::Axis)
  push!(g.axes, a)
  if length(g.axes) > prod(g.dimensions)
    g.dimensions = (length(g.axes), 1)
  end
end

function Base.push!(g::GroupPlot, p::Plot)
  push!(g, Axis(p))
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


function printObject(o::IOBuffer, object)
  print(o, "$(object)")
end

function printObject{T}(o::IOBuffer, object::AbstractArray{T,1})
  printList(o, object, brackets = true)
end


function optionHelper(o::IOBuffer, m, object; brackets=false, otherOptions=Dict{String,String}[], otherText=nothing)
  first = true
  for (sym, str) in m
    if object.(sym) != nothing
      if first
        first = false
        if brackets
          print(o, "[")
        end
      else
        print(o, ", ")
      end
      if length(str) > 0
        print(o, "$str = ")
      end
      printObject(o, object.(sym))
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

function plotHelper(o::IOBuffer, p::Histogram)
  print(o, "\\addplot+ [mark=none, $(p.style), hist={")
  optionHelper(o, histogramMap, p)
  print(o, "}] table [row sep=\\\\, y index = 0] {")
  println(o, "data\\\\")
  for d in p.data
    println(o, "$d \\\\ ")
  end
  println(o, "};")
end

function plotHelper(o::IOBuffer, p::Linear)
  print(o, "\\addplot+ ")
  optionHelper(o, linearMap, p, brackets=true)
  println(o, "coordinates {")
  for i = 1:size(p.data,2)
    println(o, "($(p.data[1,i]), $(p.data[2,i]))")
  end
  println(o, "};")
  if p.legendentry != nothing
    println(o, "\\addlegendentry{$(p.legendentry)}")
  end
end

function plotHelper(o::IOBuffer, p::Node)
  if p.style != nothing
    println(o, "\\node at (axis cs:$(p.x), $(p.y)) [$(p.style)] {$(p.data)};")
  else
    println(o, "\\node at (axis cs:$(p.x), $(p.y)) {$(p.data)};")
  end
end


function plotHelper(o::IOBuffer, p::ErrorBars)
  print(o, "\\addplot+ [")
  optionHelper(o, errorbarsMap, p)
  print(o, ",error bars/.cd, x dir=both, x explicit, y dir=both, y explicit")
  println(o,"]")
  println(o, "coordinates {")
  for i = 1:size(p.data,2)
    println(o, "($(p.data[1,i]), $(p.data[2,i])) +=($(p.data[3,i]),$(p.data[4,i])) -=($(p.data[5,i]),$(p.data[6,i]))")
  end
  println(o, "};")
  if p.legendentry != nothing
    println(o, "\\addlegendentry{$(p.legendentry)}")
  end
end

function plotHelper(o::IOBuffer, p::Quiver)
  print(o, "\\addplot+ ")
  optionHelper(o, quiverMap, p, brackets=true, otherOptions=["quiver"=>"{u=\\thisrow{u},v=\\thisrow{v}}"])
  println(o, "table {")
  println(o, "x y u v")
  for i = 1:size(p.data,2)
    println(o, "$(p.data[1,i]) $(p.data[2,i]) $(p.data[3,i]) $(p.data[4,i])")
  end
  println(o, "};")
  if p.legendentry != nothing
    println(o, "\\addlegendentry{$(p.legendentry)}")
  end
end

function plotHelper(o::IOBuffer, p::Contour)
  try
    success(`gnuplot --version`)
  catch
    error("You must have gnuplot installed and on your path to do contour plots.")
  end
  print(o, "\\addplot3[contour gnuplot={")
  optionHelper(o, contourMap, p)
  println(o, "}, mesh/cols=$(p.cols), mesh/rows=$(p.rows)] coordinates {")
  for i = 1:size(p.data,2)
    println(o, "($(p.data[1,i]), $(p.data[2,i]), $(p.data[3,i]))")
  end
  println(o, "};")
end


function plotHelper(o::IOBuffer, p::Image)
  println(o, "\\addplot [point meta min=$(p.zmin), point meta max=$(p.zmax)] graphics [xmin=$(p.xmin), xmax=$(p.xmax), ymin=$(p.ymin), ymax=$(p.ymax)] {$(p.filename)};")
end

# plot option string and contents; no \begin{axis} or \nextgroupplot
function plotHelper(o::IOBuffer, axis::Axis)
  optionHelper(o, axisMap, axis, brackets=true, otherText=[axisOptions(p) for p in axis.plots])
  for p in axis.plots
      plotHelper(o, p)
  end
end

# plot option string and contents; no \begin{axis} or \nextgroupplot
function plotHelper(o::IOBuffer, axis::PolarAxis)
  optionHelper(o, polarAxisMap, axis, brackets=true, otherText=[axisOptions(p) for p in axis.plots])
  for p in axis.plots
      plotHelper(o, p)
  end
end

function plot(axis::Axis)
  o = IOBuffer()
  print(o, "\\begin{axis}")
  plotHelper(o, axis)
  println(o, "\\end{axis}")
  TikzPicture(takebuf_string(o), options=pgfplotsoptions(), preamble=preamble)
end

function plot(axis::PolarAxis)
  o = IOBuffer()
  print(o, "\\begin{polaraxis}")
  plotHelper(o, axis)
  println(o, "\\end{polaraxis}")
  TikzPicture(takebuf_string(o), options=pgfplotsoptions(), preamble=preamble)
end

function plot(p::GroupPlot)
  o = IOBuffer()
  println(o, "\\begin{groupplot}[group style={group size=$(p.dimensions[1]) by $(p.dimensions[2])}]")
  for a in p.axes
    print(o, "\\nextgroupplot ")
    plotHelper(o, a)
  end
  println(o, "\\end{groupplot}")
  mypreamble = preamble * "\\usepgfplotslibrary{groupplots}"
  TikzPicture(takebuf_string(o), options=pgfplotsoptions(), preamble=mypreamble)
end

typealias Plottable Union(Plot,GroupPlot,Axis,PolarAxis)

plot(p::Plot) = plot(Axis(p))

plot{A<:Real,B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}) = plot(Linear(x, y))

function Plots.Linear(f::Function, range::(Real,Real); mark="none", style=nothing, legendentry=nothing)
  x = linspace(range[1], range[2])
  y = map(f, x)
  Linear(x, y, mark=mark, style=style, legendentry=legendentry)
end

plot(f::Function, range::(Real,Real)) = plot(Linear(f, range))

Base.mimewritable(::MIME"image/svg+xml", p::Plottable) = true

cleanup(p::Axis) = map(cleanup, p.plots)

cleanup(p::PolarAxis) = map(cleanup, p.plots)

cleanup(p::GroupPlot) = map(cleanup, p.axes)

cleanup(p::Plot) = nothing

cleanup(p::Image) = rm(p.filename)

axisOptions(p::Plot) = nothing

function axisOptions(p::Image)
  if p.colorbar
    return "enlargelimits = false, axis on top, colormap/blackwhite, colorbar"
  else
    return "enlargelimits = false, axis on top"
  end
end

function Base.writemime(f::IO, a::MIME"image/svg+xml", p::Plottable)
  r = Base.writemime(f, a, plot(p))
  cleanup(p)
  r
end

function save(filename::String, o::Plottable)
  _, ext = splitext(filename)
  ext = lowercase(ext)
  if ext == ".pdf"
    save(PDF(filename), plot(o))
  elseif ext == ".svg"
    save(SVG(filename), plot(o))
  elseif ext == ".tex"
    save(TEX(filename), plot(o))
  elseif ext == "." || ext == ""
    error("You must specify a file extension.")
  else
    error("Unsupported file extensions: $ext")
  end
end

end # module
