

abstract type AbstractWaveform <: Function end

(f::AbstractWaveform)(x) = evaluate(f,x)
space(::AbstractWaveform) = Fourier(PeriodicSegment(0..1))
function coefficients(f::AbstractWaveform,N::Integer)
    return [coefficient(f,n) for n in 0:N]
end

Fun(f::AbstractWaveform,N::Integer) = Fun(space(f),coefficients(f,N))
# SawtoothWave

struct SawtoothWave <: AbstractWaveform end
evaluate(::SawtoothWave,t) = t - floor(t)
function coefficient(::SawtoothWave,n::Integer)
    if n==0
        return 1/2
    elseif mod(n,2) == 1
        m = div(n+1,2)
        return -inv(m*π)
    else
        return 0
    end
end

# # SquareWave


# struct Waveform <: Function
#     f
#     coefficients
#     name::Symbol
# end
# (f::Waveform)(x) = f.f(x)
# ApproxFun.coefficients(f::Waveform,n) = f.coefficients(n)
# ApproxFun.Fun(f::Waveform,N::Integer) = Fun(Fourier(PeriodicSegment(0..1)),[coefficients(f,n) for n = 0:N])

# function sawtooth(t)
#     return t - floor(t)
# end
# function sawtoothc(n)
#     if n==0
#         return 1/2
#     elseif mod(n,2) == 1
#         m = div(n+1,2)
#         return -1/(m*pi)
#     else
#         return 0
#     end
# end

# SawtoothWave = Waveform(sawtooth,sawtoothc,:sawtooth)

# end

# # ╔═╡ 6043ac69-567c-4544-8cf7-092f232e4a20
# begin
# periodicvariation(f::Function) = Fun(z->f(z),Laurent(PeriodicSegment(0..2π)))
# periodicvariation(f::Function,m) = Fun(z->f(z),Laurent(PeriodicSegment(0..2π)),m)
# end
