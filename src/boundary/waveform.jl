import ApproxFun: Fun, coefficients, evaluate, space

abstract type AbstractWaveform <: Function end

(f::AbstractWaveform)(x) = evaluate(f,x)
space(::AbstractWaveform) = Fourier(PeriodicSegment(0..1))
coefficients(f::AbstractWaveform,N::Integer) = [coefficient(f,n) for n in 0:N]Fun(f::AbstractWaveform,N::Integer) = Fun(space(f),coefficients(f,N))
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

# SquareWave
