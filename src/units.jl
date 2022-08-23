const WavenumberDim = Unitful.𝐋^-1
const LengthDim     = Unitful.𝐋
const SpeedDim      = Unitful.𝐋*Unitful.𝐓^-1
const FrequencyDim  = Unitful.𝐓^-1
const DensityDim    = Unitful.𝐌*Unitful.𝐋^-3

# NoUnits forces simplification of μm/m = 1e-6 and Hz*s = 1

function nondimensionalize(d,x::Quantity{T,WavenumberDim,S}) where {T,S}
    return NoUnits(x*reference_length(d))
end
function nondimensionalize(d,x::Quantity{T,FrequencyDim,S}) where {T,S}
    return NoUnits(x*reference_time(d))
end
function nondimensionalize(d,x::Quantity{T,LengthDim,S}) where {T,S}
    return NoUnits(x/reference_length(d))
end
function nondimensionalize(d,x::Quantity{T,SpeedDim,S}) where {T,S}
    return NoUnits(x/reference_speed(d))
end
function nondimensionalize(d,x::Quantity{T,DensityDim,S}) where {T,S}
    return NoUnits(x/reference_density(d))
end

@with_kw struct ReferenceDimensions{T}
    length::Quantity{T,LengthDim} = 0.1u"m"
    speed::Quantity{T,SpeedDim} = 343.0u"m/s"
    density::Quantity{T,DensityDim} = 1.0u"kg*m^-3"
end

reference_length(d::ReferenceDimensions) = d.length
reference_speed(d::ReferenceDimensions) = d.speed
reference_density(d::ReferenceDimensions) = d.density
reference_time(d::ReferenceDimensions) = reference_length(d)/reference_speed(d)
