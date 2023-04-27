
abstract type Model end
masses(ℳ) = ℳ.ms # general method

struct OneMinusMinus{T} <: Model
    ms::MassTuple
    ξ::T
end
# 
function 𝒜(ℳ::OneMinusMinus, σs)
    ϕ0 = Kibble(σs, masses(ℳ)^2)
    return ϕ0 > 0 ? 0.0 : sqrt(-ϕ0) *
           (𝒜(ℳ.ξ, σs[1]) + 𝒜(ℳ.ξ, σs[2]) + 𝒜(ℳ.ξ, σs[3]))
end

struct PureSum{T} <: Model
    ms::MassTuple
    ξ::T
end

function 𝒜(ℳ::PureSum, σs)
    return 𝒜(ℳ.ξ, σs[1]) + 𝒜(ℳ.ξ, σs[2]) + 𝒜(ℳ.ξ, σs[3])
end

struct PhaseSpace <: Model
    ms::MassTuple
end
# 
function 𝒜(ℳ::PhaseSpace, σs)
    ϕ0 = Kibble(σs, masses(ℳ)^2)
    return ϕ0 > 0 ? 0.0im : 1.0 + 0.0im
end


@recipe f(ℳ::Model) = (Base.Fix1(ℐ, ℳ), masses(ℳ))


# intensity
"""
    ℐ(X, args...) = abs2(𝒜(X, args...))

A general method to calculate intensity
"""
ℐ(X, args...) = abs2(𝒜(X, args...))
ℛℯ(X, args...) = real(𝒜(X, args...))
ℐ𝓂(X, args...) = imag(𝒜(X, args...))
