# Mixing length model

# Equation 2.38
function mixing_length(δ,x)
    Δ = 0.085*δ
    return Δ*tanh(0.41*x/Δ)
end

# Equation 2.39
mixing_length2integral_length(l) = l/0.41
integral_length2mixing_length(Λ) = 0.41*Λ


# Source term

# Equation 2.17

k = sqrt(k₁^2 + k₃^2)

integrand(x₂,y₂) = (k₁/k)^2 *exp(-(x₂+y₂)*k)*dU₁(x₂)*dU₁(y₂)*φ₂₂(x₂,y₂)

bertagnolio2014_stretch(γ) = @SVector [ 0.4; γ^(1/5); sqrt(2γ)]
stalnov2016_stretch = @SVector [1.0; 0.5; 0.75]
fischer2017_stretch = @SVector [1.0; 0.74; 0.9]
