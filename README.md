# PGFPlots

This library uses the LaTeX package [pgfplots](http://ctan.org/pkg/pgfplots) to produce plots. It integrates with IJulia, outputting SVG images to the notebook.

In order to use this package, you need to install TikzPictures by running `Pkg.clone("https://github.com/sisl/TikzPictures.jl)`. You also need to have the pgfplots package (version 1.10 or later) installed.

## Examples

### Basic plot
```julia
using PGFPlots
a = [1, 2, 3]
b = [1, 3, 2]
plot(a, b)
```

### Group plots
```julia
using PGFPlots
using Distributions
# define a normal distribution with mean 3 and standard deviation 2
dist = Normal(3,2)
# generate 1000 samples
d = rand(dist, 1000)
# try out histograms with a variety of different number of bins
bins = [3 10 50 1000]
g = GroupPlot(2,2)
for i = 1:length(bins)
    push!(g, Axis(Plots.Histogram(d, bins=bins[i]), ymin=0))
end
g
```

### Function plots

```julia
using PGFPlots
plot(sin,(0,10))
```
