using PGFPlots
using Base.Test

@assert success(`lualatex -v`)
Pkg.add("NBInclude")
using NBInclude
nbinclude(Pkg.dir("PGFPlots", "doc", "PGFPlots.ipynb"))
