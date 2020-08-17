using DataFrames

Plots.Linear(df::DataFrame; x::Symbol=:x, y::Symbol=:y, kwargs...) =
    Plots.Linear(df[!,x], df[!,y]; kwargs...)

Plots.Linear3(df::DataFrame; x::Symbol=:x, y::Symbol=:y, z::Symbol=:z, kwargs...) =
    Plots.Linear3(df[!,x], df[!,y], df[!,z]; kwargs...)

Plots.Histogram(df::DataFrame; x::Symbol=:x, kwargs...) =
    Plots.Histogram(df[!,x]; kwargs...)

Plots.Scatter(df::DataFrame; x::Symbol=:x, y::Symbol=:y, kwargs...) =
    Plots.Scatter(df[!,x], df[!,y]; kwargs...)

Plots.Quiver(df::DataFrame; x::Symbol=:x, y::Symbol=:y, u::Symbol=:u, v::Symbol=:v, kwargs...) =
    Plots.Quiver(convert(Vector, df[!,x]), convert(Vector, df[!,y]), convert(Vector, df[!,u]), convert(Vector, df[!,v]); kwargs...)

Plots.Histogram2(df::DataFrame; x::Symbol=:x, y::Symbol=:y, kwargs...) =
    Plots.Histogram2(convert(Vector, df[!,x]), convert(Vector, df[!,y]); kwargs...)

Plots.Histogram2(df::DataFrame, edges_x::AbstractVector{C}, edges_y::AbstractVector{C}; x::Symbol=:x, y::Symbol=:y, kwargs...) where {C<:Real} =
    Plots.Histogram2(convert(Vector, df[!,x]), convert(Vector, df[!,y]), edges_x, edges_y; kwargs...)
