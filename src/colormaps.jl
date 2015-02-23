module ColorMaps

export ColorMap, Gray, RGBArray, Distinguishable, write
import Images: grayim, imwrite, ImageCmap
import Color: RGB, distinguishable_colors

abstract ColorMap

type Gray <: ColorMap
end

type RGBArray <: ColorMap
    colors::Array{RGB{Float64},1}
end

Distinguishable(n::Integer) = RGBArray(convert(Array{RGB{Float64},1}, distinguishable_colors(n)))

function Base.write(colormap::ColorMap, data, filename)
    imwrite(grayim(data), filename)
end

function Base.write(colormap::ColorMaps.RGBArray, data, filename)
    img = ImageCmap(uint8(1.+(length(colormap.colors)-1).*data), colormap.colors)
    imwrite(img, filename)
end

end
