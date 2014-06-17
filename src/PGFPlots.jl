module PGFPlots

export plot, Axis

using TikzPictures

type Plot
  data
end

type Axis
  plots
  title::String
  xlabel::String
  ylabel::String
  Axis(data; title="", xlabel="", ylabel="") = new([Plot(data)],
  title, xlabel, ylabel
  )
end

function optionString(axis::Axis)
  t = ""
  first = true
  for n in names(Axis)
    if n != :plots && length(axis.(n)) > 0
      if first
        first = false
      else
        t *= ", "
      end
      t *= "$(string(n)) = $(axis.(n))"
    end
  end
  if length(t) > 0
    t = "[$t]"
  else
    return ""
  end
end

function plot(axis::Axis)
  t = "\\begin{axis}$(optionString(axis))\n"
  for p in axis.plots
    t *= "\\addplot coordinates {\n"
    for i = 1:size(p.data,2)
      t *= "($(p.data[1,i]), $(p.data[2,i]))\n"
    end
    t *= "};\n"
  end
  t *= "\\end{axis}"
  preamble = "\\usepackage{pgfplots}\n"
  preamble *= "\\pgfplotsset{compat=newest}"
  TikzPicture(t, options="scale=1", preamble=preamble)
end

plot(data) = plot(Axis(data))

plot(x::AbstractArray, y::AbstractArray) = plot([x y]')


end # module
