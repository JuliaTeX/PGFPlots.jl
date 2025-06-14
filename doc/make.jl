using Documenter, PGFPlots, Reel

makedocs(sitename="PGFPlots Docs", remotes = nothing, format=Documenter.HTML())

deploydocs(
    repo = "github.com/JuliaTeX/PGFPlots.jl.git", 
    push_preview=true
)