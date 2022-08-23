#module PlaneWaves



"""
    wave = PlaneWave{N,T} <: Function

    represents plane wave

        x -> Aexp(ikx - iωt)

"""

struct PlaneWave{N,T} <: Function
    k::SVector{N,T}
    ω::T
    A::T
    function PlaneWave(k::SVector{N,S}, ω=0, A=1) where {N,S}
        T = eltype(promote(k[1], ω, A))
        new{N,T}(convert.(T, k), convert(T, ω), convert(T, A))
    end
end

wavevector(a::PlaneWave) = a.k
frequency(a::PlaneWave) = a.ω
amplitude(a::PlaneWave) = a.A

wavenumber(a::PlaneWave, i) = getindex(wavevector(a), i)

setamplitude(a::PlaneWave, λ::Number) = PlaneWave(wavevector(a), frequency(a), λ)

/(a::PlaneWave, b::Number) = setamplitude(a, amplitude(a) / b)
*(a::PlaneWave, b::Number) = setamplitude(a, amplitude(a) * b)
-(a::PlaneWave) = setamplitude(a, -amplitude(a))
*(a::Number, b::PlaneWave) = *(b, a)

spacederivative(a::PlaneWave, n::Integer, order::Integer) = a * (im * wavenumber(a, n))^order

function evaluate(wave::PlaneWave, x::SVector, t=0)
    @unpack k, ω, A = wave
    return A * exp(im * sum(k .* x) - im * ω * t)
end
function evaluate(wave::PlaneWave, z::Complex, t=0)
    @unpack k, ω, A = wave
    x = SVector(reim(z))
    return A * exp(im * sum(k .* x) - im * ω * t)
end

(wave::PlaneWave)(args...) = evaluate(wave, args...)

islefthalfline(a,b) = (a == -Inf) && isfinite(b)
isrighthalfline(a,b) = isfinite(a) && (b == Inf)
isfiniteinterval(a,b) =  isfinite(a) && isfinite(b)

function fourierx(wave::PlaneWave,a::Real,b::Real,shift=0)
    A = amplitude(wave);
    δ = wavenumber(wave,1);
    r = ExponentialShift(ScalarPole(-δ,-A/im),b-shift,-b*δ)
    l = ExponentialShift(ScalarPole(-δ, A/im),a-shift,-a*δ)
    if islefthalfline(a,b)  return r end
    if isrighthalfline(a,b) return l end
    return r-l
end


function fourierx(wave::PlaneWave,I::Interval,shift=0)
    fourierx(wave,I.left,I.right,shift)
end


#end
