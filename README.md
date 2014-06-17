# PGFPlots

This library uses the LaTeX package [pgfplots](http://ctan.org/pkg/pgfplots) to produce plots. It integrates with IJulia, outputting SVG images to the notebook.

In order to use this package, you need to install TikzPictures by running `Pkg.clone("https://github.com/sisl/TikzPictures.jl)`. You also need to have the pgfplots package installed.

## Example

```julia
using PGFPlots
a = [1, 2, 3]
b = [1, 3, 2]
plot(a, b)
```
