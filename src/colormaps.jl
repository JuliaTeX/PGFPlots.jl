module ColorMaps

export ColorMap, Gray, Brew, RGBArray, Distinguishable, SparseDistinguishable, write
import Images: grayim, save, ImageCmap
import Colors: RGB, distinguishable_colors, colormap
import ColorBrewer: palette

abstract ColorMap

type Gray <: ColorMap
    invert::Bool
    Gray(;invert = false) = new(invert)
end

type RGBArray <: ColorMap
    colors::Array{RGB{Float64},1}
    function RGBArray(colors; invert=false)
        if invert
            colors = flipdim(colors, 1)
        end
        new(colors)
    end
end

function Brew(name::AbstractString, number::Integer; invert = false)
    a = RGB{Float64}[convert(RGB{Float64}, c) for c in palette(name, number)]
    RGBArray(a, invert=invert)
end

function Distinguishable(n::Integer; invert=false)
    RGBArray(convert(Array{RGB{Float64},1}, distinguishable_colors(n)), invert=invert)
end

function SparseDistinguishable(n, nonzeroidx; zerocolor = RGB{Float64}(1.,1.,1.))
    m = RGB{Float64}[zerocolor for i = 1:n]
    sm = ColorMaps.Distinguishable(length(nonzeroidx))
    for i = 1:length(sm.colors)
        m[nonzeroidx[i]] = sm.colors[i]
    end
    ColorMaps.RGBArray(m)
end

function Named(name::AbstractString = "Jet", levels=255)
    if isequal(name, "Jet")
        return RGBArray([RGB(min(max(min(4*i/256-1.5,-4*i/256+4.5),0),1), min(max(min(4*i/256-0.5,-4*i/256+3.5),0),1),min(max(min(4*i/256+0.5,-4*i/256+2.5),0),1)) for i = 1:255])
    else
        return RGBArray(colormap(name, levels))
    end
end

function Base.write(colormap::ColorMaps.Gray, data, filename)
    if colormap.invert
        inverteddata = 1. - data
        save(filename, grayim(rotr90(flipdim(inverteddata, 1))))
    else
        save(filename, grayim(rotr90(flipdim(data, 1))))
    end
end

Base.write(colormap::ColorMap, data, filename) = error("Not supported")

function Base.write(colormap::ColorMaps.RGBArray, data, filename)
    img = ImageCmap(round(UInt8, 1.+(length(colormap.colors)-1).*(data)), colormap.colors)
    save(filename, img)
end

end
