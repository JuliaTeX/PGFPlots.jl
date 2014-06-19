module PGFPlots

export plot, Axis, GroupPlot, Plots

include("ndgrid.jl")

import Images: grayim, imwrite

preamble = """
\\usepackage{pgfplots}
\\pgfplotsset{compat=newest}
"""

using TikzPictures

module Plots

export Plot, Histogram, Linear, Image

abstract Plot

type Histogram <: Plot
  data::AbstractArray{Real,1}
  bins::Integer
  density::Bool
  cumulative::Bool
  style::String
  Histogram(data; bins=20, density=false, cumulative=false, style="fill=blue!10") = new(data, bins, density, cumulative, style)
end

type Linear <: Plot
  data::AbstractArray{Real,2}
  mark
  style
  Linear(data::AbstractArray{Real,2}; mark=nothing, style=nothing) = new(data, mark, style)
  Linear{A<:Real, B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}; mark=nothing, style=nothing) = new([x y]', mark, style)
end

type Image <: Plot
  filename::String
  xmin::Real
  xmax::Real
  ymin::Real
  ymax::Real
end

end # end plot module

histogramMap = [
  :bins => "bins",
  :density => "density",
  :cumulative => "cumulative"
  ]

linearMap = [
  :mark => "mark",
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
  Axis(plot::Plot;title=nothing, xlabel=nothing, ylabel=nothing, xmin=nothing, xmax=nothing,
       ymin=nothing, ymax=nothing, enlargelimits=nothing, axisOnTop=nothing) =
    new([plot], title, xlabel, ylabel, xmin, xmax, ymin, ymax, enlargelimits, axisOnTop
  )
end

axisMap = [
  :title => "title",
  :xlabel => "xlabel",
  :ylabel => "ylabel",
  :xmin => "xmin",
  :ymin => "ymin",
  :enlargelimits => "enlargelimits",
  :axisOnTop => "axis on top"
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
        print(o, "$str = $(object.(sym))")
      else
        print(o, "$(object.(sym))")
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

function Plots.Image(f::Function, xrange::(Real,Real), yrange::(Real,Real); filename="tmp.png")
  x = linspace(xrange[1], xrange[2])
  y = linspace(yrange[1], yrange[2])
  (X, Y) = meshgrid(x, y)
  A = map(f, X, Y)
  # normalize A
  A = A .- minimum(A)
  A = A ./ maximum(A)
  A = rotr90(A)
  imwrite(grayim(A), filename)
  Image(filename, xrange[1], xrange[2], yrange[1], yrange[2])
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
