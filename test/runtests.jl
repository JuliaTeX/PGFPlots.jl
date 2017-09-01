using PGFPlots
using Base.Test

using Reel
Reel.extension(m::MIME"image/svg+xml") = "svg"

@assert success(`lualatex -v`)
using NBInclude
nbinclude(joinpath(dirname(@__FILE__), "..", "doc", "PGFPlots.ipynb"))
