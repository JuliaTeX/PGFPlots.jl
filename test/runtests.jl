using PGFPlots
using Base.Test

@assert success(`lualatex -v`)
using NBInclude
nbinclude(joinpath(dirname(@__FILE__), "..", "doc", "PGFPlots.ipynb"))
