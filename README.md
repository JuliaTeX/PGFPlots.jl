# PGFPlots

This library uses the LaTeX package [pgfplots](http://ctan.org/pkg/pgfplots) to produce plots. It integrates with IJulia, outputting SVG images to the notebook.

In order to use this package, you need the following:
* TikzPictures. Obtain by running `Pkg.clone("https://github.com/sisl/TikzPictures.jl")`. TikzPictures requires pdf2svg. On Ubuntu, you can get this by running `sudo apt-get install pdf2svg`. On Windows, you can download the binaries from http://www.cityinthesky.co.uk/opensource/pdf2svg/. Be sure to add pdf2svg to your path (and restart).
* Pgfplots (version 1.10 or later). Install using your latex package manager.

## Examples

### Basic plot
```julia
using PGFPlots
a = [1, 2, 3]
b = [1, 3, 2]
plot(a, b)
```

### Group plots (with histograms)
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

### 2D image plots

```julia
using PGFPlots
Plots.Image((x,y)->sin(x)*y, (0,10), (0,10))
```
