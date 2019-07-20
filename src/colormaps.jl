module ColorMaps

using ColorSchemes
export ColorMap, GrayMap, Brew, RGBArrayMap, Distinguishable, SparseDistinguishable, write
import Images: colorview, save, Gray, ImageMeta, clamp01nan 
import Colors: RGB, distinguishable_colors, colormap
import ColorBrewer: palette
import IndirectArrays: IndirectArray

abstract type ColorMap end

mutable struct GrayMap <: ColorMap
    invert::Bool
    GrayMap(;invert = false) = new(invert)
end

mutable struct RGBArrayMap <: ColorMap
    colors::Vector{RGB{Float64}}
    interpolation_levels::UInt
    function RGBArrayMap(colors::AbstractVector{RGB{Float64}}; invert=false, interpolation_levels=0)
        if invert
            colors = reverse(colors)
        end
        new(colors, interpolation_levels)
    end
    function RGBArrayMap(colors::ColorScheme; invert=false, interpolation_levels=0)
        return RGBArrayMap(colors.colors, invert=invert, interpolation_levels=interpolation_levels)
    end
end

function Brew(name::AbstractString, number::Integer; invert = false)
    a = RGB{Float64}[convert(RGB{Float64}, c) for c in palette(name, number)]
    RGBArrayMap(a, invert=invert)
end

function Distinguishable(n::Integer; invert=false)
    RGBArrayMap(convert(Array{RGB{Float64},1}, distinguishable_colors(n)), invert=invert)
end

function SparseDistinguishable(n, nonzeroidx; zerocolor = RGB{Float64}(1.,1.,1.))
    m = RGB{Float64}[zerocolor for i = 1:n]
    sm = ColorMaps.Distinguishable(length(nonzeroidx))
    for i = 1:length(sm.colors)
        m[nonzeroidx[i]] = sm.colors[i]
    end
    ColorMaps.RGBArrayMap(m)
end

function Named(name::AbstractString = "Jet", levels=255)
    if isequal(name, "Jet")
        return RGBArrayMap([RGB(min(max(min(4*i/256-1.5,-4*i/256+4.5),0),1), min(max(min(4*i/256-0.5,-4*i/256+3.5),0),1),min(max(min(4*i/256+0.5,-4*i/256+2.5),0),1)) for i = 1:255])
    else
        return RGBArrayMap(colormap(name, levels))
    end
end

function Base.write(colormap::GrayMap, data, filename)
    if colormap.invert
        inverteddata = 1 .- data
        save(filename, colorview(Gray, inverteddata))
    else
        save(filename, colorview(Gray, data))
    end
end

Base.write(colormap::ColorMap, data, filename) = error("Not supported")

function interpolate_RGBArrayMap(colormap::RGBArrayMap)
    levels = colormap.interpolation_levels
    if levels == 0
        return colormap.colors
    end

    n = length(colormap.colors)
    colors = Vector{RGB{Float64}}(undef,levels)
    for (i,x) in enumerate(range(0.0,stop=1.0,length=levels))
        t = 1 + x*(n-1)
        a = colormap.colors[floor(Int, t)]
        b = colormap.colors[ceil(Int, t)]
        α = ceil(Int, t) - t
        colors[i] = α*a + (1-α)*b
    end
    return colors
end

function Base.write(colormap::RGBArrayMap, data, filename)
    colors = interpolate_RGBArrayMap(colormap)
    n = length(colors)
    data = clamp01nan.(data)
    if n <= 2^8
        img = ImageMeta(IndirectArray([round(UInt8, (n-1)*v + 1) for v in data], colors))
        save(filename, img)
    elseif n <= 2^16
        img = ImageMeta(IndirectArray([round(UInt16, (n-1)*v + 1) for v in data], colors))
        save(filename, img)
    elseif n <= 2^32
        img = ImageMeta(IndirectArray([round(UInt32, (n-1)*v + 1) for v in data], colors))
        save(filename, img)
    else
        error("ColorMap has too many colors")
    end
end

end
