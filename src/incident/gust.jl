function convectedgust(;M=0.1,ω=1,k₂=1,k₃=0,A=1)
	β = sqrt(1-M^2)
	k₁ = ω/(β*M)
    k = @SVector [k₁,k₂,k₃]
    return PlaneWave(k,ω,A)
end
