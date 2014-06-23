module PGFPlots

export plot, Axis, GroupPlot, Plots

include("plots.jl")

import TikzPictures: TikzPicture, PDF, TEX, SVG, save

preamble = """
\\usepackage{pgfplots}
\\pgfplotsset{compat=newest}
"""

histogramMap = [
  :bins => "bins",
  :density => "density",
  :cumulative => "cumulative"
  ]

linearMap = [
  :mark => "mark",
  :style => ""
  ]

contourMap = [
  :number => "number",
  :levels => "levels",
  :style => ""
]


using .Plots

type Axis
  plots::AbstractArray{Plot,1}
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

  Axis(plot::Plot;title=nothing, xlabel=nothing, ylabel=nothing, xmin=nothing, xmax=nothing,
       ymin=nothing, ymax=nothing, enlargelimits=nothing, axisOnTop=nothing, view="{0}{90}") =
    new([plot], title, xlabel, ylabel, xmin, xmax, ymin, ymax, enlargelimits, axisOnTop, view
  )
end

axisMap = [
  :title => "title",
  :xlabel => "xlabel",
  :ylabel => "ylabel",
  :xmin => "xmin",
  :ymin => "ymin",
  :enlargelimits => "enlargelimits",
  :axisOnTop => "axis on top",
  :view => "view"
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


function optionHelper(o::IOBuffer, m, object; brackets=false)
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
    print(o, "$d \\\\ ")
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
  println(o, "\\addplot graphics [xmin=$(p.xmin), xmax=$(p.xmax), ymin=$(p.ymin), ymax=$(p.ymax)] {$(p.filename)};")
end

# plot option string and contents; no \begin{axis} or \nextgroupplot
function plotHelper(o::IOBuffer, axis::Axis)
  optionHelper(o, axisMap, axis, brackets=true)
  for p in axis.plots
      plotHelper(o, p)
  end
end

function plot(axis::Axis)
  o = IOBuffer()
  print(o, "\\begin{axis}")
  plotHelper(o, axis)
  println(o, "\\end{axis}")
  TikzPicture(takebuf_string(o), options="scale=1", preamble=preamble)
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
  TikzPicture(takebuf_string(o), options="scale=1", preamble=mypreamble)
end

plot(p::Plot) = plot(Axis(p))

plot(p::Image) = plot(Axis(p, enlargelimits=false, axisOnTop=true))

plot{A<:Real,B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}) = plot(Linear(x, y))

function Plots.Linear(f::Function, range::(Real,Real))
  x = linspace(range[1], range[2])
  y = map(f, x)
  Linear(x, y, mark="none")
end

plot(f::Function, range::(Real,Real)) = plot(Linear(f, range))

Base.mimewritable(::MIME"image/svg+xml", p::Plot) = true
Base.writemime(f::IO, a::MIME"image/svg+xml", p::Plot) = Base.writemime(f, a, plot(p))

Base.mimewritable(::MIME"image/svg+xml", p::GroupPlot) = true
Base.writemime(f::IO, a::MIME"image/svg+xml", p::GroupPlot) = Base.writemime(f, a, plot(p))

Base.mimewritable(::MIME"image/svg+xml", p::Axis) = true
Base.writemime(f::IO, a::MIME"image/svg+xml", p::Axis) = Base.writemime(f, a, plot(p))


end # module
