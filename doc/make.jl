using Documenter, PGFPlots, Reel

makedocs(sitename="PGFPlots Docs", remotes = nothing, format=Documenter.HTML(), pages=["Home" => "index.md", "Installation" => "installation.md", "Examples" => "examples.md", "Future Plans" => "future_plans.md"])

deploydocs(
    repo = "github.com/JuliaTeX/PGFPlots.jl.git", 
    push_preview=true
)