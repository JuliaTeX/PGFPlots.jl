module Plots

export Plot, Histogram, Linear, Image, Contour, Scatter
import Images: grayim, imwrite

include("ndgrid.jl")

abstract Plot

type Histogram <: Plot
  data::AbstractArray{Real,1}
  bins::Integer
  density::Bool
  cumulative::Bool
  style::String
  Histogram(data; bins=20, density=false, cumulative=false, style="fill=blue!10") = new(data, bins, density, cumulative, style)
end

type Contour <: Plot
  data::AbstractArray{Real,2} # 3 x n matrix
  cols::Integer
  rows::Integer
  style
  number
  levels
  Contour(data, cols, rows; style=nothing, number=nothing, levels=nothing) = new(data, cols, rows, style, number, levels)
  function Contour(f::Function, xrange::(Real,Real), yrange::(Real,Real); style=nothing, number=nothing, levels=nothing)
    x = linspace(xrange[1], xrange[2], 40)
    y = linspace(yrange[1], yrange[2], 40)
    (X, Y) = meshgrid(x, y)
    A = map(f, X, Y)
    A = [X[:]'; Y[:]'; A[:]']
    new(A, length(x), length(y), style, number, levels)
  end
end

type Linear <: Plot
  data::AbstractArray{Real,2}
  mark
  style
  legendentry
  onlyMarks
  Linear(data::AbstractArray{Real,2}; mark=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = new(data, mark, style, legendentry, onlyMarks)
end

Linear{A<:Real, B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}; mark=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = Linear([x y]', mark=mark, style=style, legendentry=legendentry, onlyMarks=onlyMarks)
Linear{A<:Real}(data::AbstractArray{A,1}; mark=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = Linear([1:length(data)], data, mark=mark, style=style, legendentry=legendentry, onlyMarks=onlyMarks)
Linear{T<:Real}(data::AbstractArray{T,2}; mark=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = Linear(data, mark=mark, style=style, legendentry=legendentry, onlyMarks=onlyMarks)
Linear{A<:Real, B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}; mark=nothing, style=nothing, legendentry=nothing, onlyMarks=nothing) = Linear([x y]', mark=mark, style=style, legendentry=legendentry, onlyMarks=onlyMarks)

Scatter{T<:Real}(data::AbstractArray{T,2}; mark=nothing, style=nothing, legendentry=nothing) = Linear(data, mark=mark, style=style, onlyMarks = true, legendentry=nothing)
Scatter{A<:Real, B<:Real}(x::AbstractArray{A,1}, y::AbstractArray{B,1}; mark=nothing, style=nothing, legendentry=nothing) = Scatter([x y]', mark=mark, style=style, legendentry=nothing)
Scatter{A<:Real, B<:Real}(x::A, y::B; mark=nothing, style=nothing, legendentry=nothing) = Scatter([x y]', mark=mark, style=style, legendentry=nothing)

global _imgid = 1

type Image <: Plot
  filename::String
  xmin::Real
  xmax::Real
  ymin::Real
  ymax::Real
  function Image(A::Matrix{Float64}, xrange::(Real,Real), yrange::(Real,Real); filename=nothing)
    global _imgid
    if filename == nothing
      filename = "tmp_$(_imgid).png"
      _imgid += 1
    end
    A = A .- minimum(A)
    A = A ./ maximum(A)
    imwrite(grayim(A), filename)
    new(filename, xrange[1], xrange[2], yrange[1], yrange[2])
  end
  function Image(f::Function, xrange::(Real,Real), yrange::(Real,Real); filename=nothing)
    x = linspace(xrange[1], xrange[2])
    y = linspace(yrange[1], yrange[2])
    (X, Y) = meshgrid(x, y)
    A = map(f, X, Y)
    A = rotr90(A)
    Image(A, xrange[1], xrange[2], yrange[1], yrange[2], filename=filename)
  end

end

end # end plot module
