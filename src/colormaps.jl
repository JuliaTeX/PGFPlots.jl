module ColorMaps

export ColorMap, Gray, RGBArray, Distinguishable, SparseDistinguishable, write
import Images: grayim, imwrite, ImageCmap
import Colors: RGB, distinguishable_colors, colormap

abstract ColorMap

type Gray <: ColorMap
    invert::Bool
    Gray(;invert = false) = new(invert)
end

type RGBArray <: ColorMap
    colors::Array{RGB{Float64},1}
end

Distinguishable(n::Integer) = RGBArray(convert(Array{RGB{Float64},1}, distinguishable_colors(n)))

function SparseDistinguishable(n, nonzeroidx; zerocolor = RGB{Float64}(1.,1.,1.))
    m = RGB{Float64}[zerocolor for i = 1:n]
    sm = ColorMaps.Distinguishable(length(nonzeroidx))
    for i = 1:length(sm.colors)
        m[nonzeroidx[i]] = sm.colors[i]
    end
    ColorMaps.RGBArray(m)
end

function Named(name::String = "Jet", levels=255)
    if isequal(name, "Jet")
        return RGBArray([RGB(min(max(min(4*i/256-1.5,-4*i/256+4.5),0),1), min(max(min(4*i/256-0.5,-4*i/256+3.5),0),1),min(max(min(4*i/256+0.5,-4*i/256+2.5),0),1)) for i = 1:255])
    else
        return RGBArray(colormap(name, levels))
    end
end

function Base.write(colormap::ColorMaps.Gray, data, filename)
    if colormap.invert
        inverteddata = 1. - data
        imwrite(grayim(rotr90(flipud(inverteddata))), filename)
    else
        imwrite(grayim(rotr90(flipud(data))), filename)
    end
end

Base.write(colormap::ColorMap, data, filename) = error("Not supported")

function Base.write(colormap::ColorMaps.RGBArray, data, filename)
    img = ImageCmap(uint8(1.+(length(colormap.colors)-1).*(data)), colormap.colors)
    imwrite(img, filename)
end

end
