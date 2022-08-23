
export partitionline

abstract type AbstractBoundaryCondition end

struct NeumannCondition <: AbstractBoundaryCondition end
struct DirichletCondition <: AbstractBoundaryCondition end
@with_kw struct RobinCondition{T} <: AbstractBoundaryCondition
    a::T
    b::T = one(T)
end

# struct LinearCondition{T} <: AbstractBoundaryCondition
#     a::Polynomial{T}
#     b::T
# end

# A boundary is formed of a sequence of domains and boundary conditions
struct DomainCondition{T,S}
    domain::T
    boundarycondition::S
end
domain(x::DomainCondition) = x.domain
boundarycondition(x::DomainCondition) = x.boundarycondition

neumanncondition(d::Domain) = DomainCondition(d, NeumannCondition())
dirichletcondition(d::Domain) = DomainCondition(d, DirichletCondition())
robincondition(d::Domain, a, b) = DomainCondition(d, RobinCondition(a, b))

function dirichletneumann(d::Domain, dirichletfirst::Bool)
    bc(i, d) = xor(isodd(i),!dirichletfirst) ? dirichletcondition(d) : neumanncondition(d)
    vcat(bc(i,d) for (i,d) in enumerate(d))
end
function dirichletneumann(x::AbstractVector, dirichletfirst::Bool)
    dirichletneumann(partitionline(x), dirichletfirst)
end

function isincreasing(x)
    for i in eachindex(x)
        if issubset((i, i + 1), eachindex(x))
            x[i] > x[i+1] && return false
        end
    end
    return true
end

function partitionline(x::AbstractVector{<:AbstractFloat})
    @assert isincreasing(x)
    f = (first(x) !== -Inf)
    l = (last(x) !== Inf)
    v = similar(x, length(x) + f + l)
    v[begin] = -Inf
    v[begin+1:end-1] .= x[(begin+1-f):(end-1+l)]
    v[end] = Inf
    return UnionDomain(v[i] .. v[i+1] for i = firstindex(v):lastindex(v)-1)
end

function domain(a::Vector{<:DomainCondition})
    UnionDomain(domain(d) for d in a)
end

function boundarycondition(a::Vector{<:DomainCondition})
    [boundarycondition(d) for d in a]
end

# Convenience

function rigidleadingedge(x=zeros(Float64, 1))
    @assert length(x) == 1
    dirichletneumann(x, true)
end

function rigidtrailingedge(x=zeros(Float64, 1))
    @assert length(x) == 1
    dirichletneumann(x, false)
end

function rigidfiniteplate(x=[0.0, 1.0])
    @assert length(x) == 2
    dirichletneumann(x, true)
end

function slottedleadingedge(x=[0.0, 1.0, 2.0])
    @assert length(x) == 3
    dirichletneumann(x, true)
end
