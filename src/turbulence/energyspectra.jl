# Abstract energy spectrum type

abstract type AbstractEnergySpectrum{homogeneous} <: Function end

characteristic_length(Π::AbstractEnergySpectrum{true}) = Π.l
characteristic_length(Π::AbstractEnergySpectrum{true},x) = Π.l
characteristic_length(Π::AbstractEnergySpectrum{false},x) = Π.l(x)
characteristic_length(Π::AbstractEnergySpectrum,x,y) = characteristic_length(Π,geometricmean(x,y))

function xspectral_density(Π::AbstractEnergySpectrum,u,wavevector,x,y)
    u(geometricmean(x,y))*nondim_xspectral_density(Π,wavevector,x,y)
end
function autospectrum(Π::AbstractEnergySpectrum,u,wavevector,x)
    u(geometricmean(x,y))*nondim_autospectrum(Π,wavevector,x)
end

# GaussianEnergySpectrum

struct GaussianEnergySpectrum{homogeneous,T} <: AbstractEnergySpectrum{homogeneous}
    l::T
end
# Equation 2.21
function gaussian_longitudinal_correlation(l,r)
    exp(-(r/l)^2)
end
# Equation 2.22
function gaussian_longitudinal_integral_lengthscale(l)
    return sqrt(π)*l/2
end
function longitudinal_integral_lengthscale(Π::GaussianEnergySpectrum...)
    l = characteristic_length(Π...)
    gaussian_longitudinal_integral_lengthscale(l)
end
function nondim_xspectral_density(wavevector,Π::GaussianEnergySpectrum,x,y)
    l = characteristic_length(Π,x,y)
    k² = LinearAlgebra.norm(wavevector,2)^2
    return l^4*k²*exp(k²*l^2/4 - (x-y)^2/l^2)/(16π)
end

# TODO:gaussian_autocorrelation

# Vonkarman

struct VonKarmanEnergySpectrum{homogeneous,T} <: AbstractEnergySpectrum{homogeneous}
    ν::T
    l::T
end

vonkarman_parameter(Π::AbstractEnergySpectrum{true}) = Π.ν
vonkarman_parameter(Π::AbstractEnergySpectrum{true},x) = Π.ν
vonkarman_parameter(Π::AbstractEnergySpectrum{false},x) = Π.ν(x)

vonkarman_parameter(Π::VonKarmanEnergySpectrum,x,y) = vonkarman_parameter(Π,geometricmean(x,y))

function VonKarmanEnergySpectrum(homogeneous::Bool,ν,l)
    x = promote(ν,l)
    VonKarmanEnergySpectrum{homogeneous,eltype(x)}(x...)
end

# Gaussian ν = 1/3
# Liepmann ν = 1/2
# RDT ν = 7/6

# Equation 2.23
function vonkarman_longitudinal_correlation(ν,l,r)
    inv(gamma(ν)*2^(ν-1))*(r/l)^ν*besselk(ν,r/l)
end

function longitudinal_correlation(Π::VonKarmanEnergySpectrum,r)
    @unpack ν, l = Π
    inv(gamma(ν)*2^(ν-1))*(r/l)^ν*besselk(ν,r/l)
end

# Equation 2.24
function vonkarman_longitudinal_integral_lengthscale(ν,l)
    return (sqrt(π)*SpecialFunctions.gamma(ν+1/2)*l)/SpecialFunctions.gamma(ν)
end
function longitudinal_integral_lengthscale(Π::VonKarmanEnergySpectrum...)
    ν = vonkarman_parameter(Π...)
    l = characteristic_length(Π...)
    vonkarman_longitudinal_integral_lengthscale(ν,l)
end

# Equation 2.26
normalised_zeta(l,k,x,y) = sqrt(1+(l*k)^2)*LinearAlgebra.norm(x-y,2)/l

function nondim_xspectral_density(Π::VonKarmanEnergySpectrum,wavevector,x,y)
    l = characteristic_length(Π,x,y)
    ν = vonkarman_parameter(Π,x,y)
    k = LinearAlgebra.norm(wavevector,2)
    ζ = normalised_zeta(l,k,x,y)
    κ = l*k
    if iszero(ζ)
        ζ = 1e-12
    end
    den = (SpecialFunctions.gamma(ν)*π*(2^(ν+1))*(1+κ^2)^(ν+2))
    num = (l^2)*κ^2*ζ^(ν+2)*SpecialFunctions.besselk(ν+2,ζ)
    return num/den
end

# Equation 2.27
function nondim_autospectrum(Π,wavevector,x)
    l = characteristic_length(Π,x)
    ν = vonkarman_parameter(Π,x)
    k = LinearAlgebra.norm(wavevector,2)
    num = ν*(ν+1)*l^4*k^2
    den = (π*(1+(k*l)^2))^(ν+2)

    # small k: k^2
    # large k: k^(-2ν-2)
    # and k looks like sqrt((ω/U)^2 + k₃^2)
    return num/den
end

# Stretched spectra

struct StretchedSpectrum{homogeneous,T,S} <: AbstractEnergySpectrum{homogeneous}
    Π::T
    β::S
    function StretchedSpectrum(Π::AbstractEnergySpectrum{homogeneous},β::S) where {homogeneous,S}
        new{homogeneous,typeof(Π),S}(Π,β)
    end
end

isotropic(Π::StretchedSpectrum) = Π.Π
isotropic(Π::StretchedSpectrum,args...) = (isotropic(Π),args...)

# Equation 2.33
function nondim_autospectrum(wavevector,Π::StretchedSpectrum...)
    @unpack β = Π
    β[1]*β[3]*nondim_autospectrum(β.*wavevector,isotropic(Π...)...)
end

# Equation 2.34
function nondim_vertical_xspectral_density(wavevector,Π::StretchedSpectrum...)
    @unpack β = Π
    β[1]*β[2]*β[3]*nondim_vertical_xspectral_density(β.*wavevector,isotropic(Π...)...)
end
