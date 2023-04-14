
# S-wave BW amplitude
struct BW
    m::Float64
    Γ::Float64
end
𝒜(ξ::BW, σ::Real) = ξ.m*ξ.Γ/(ξ.m^2-σ-1im*ξ.m*ξ.Γ)
# 



struct LookupIsobar
    itr::AbstractInterpolation
end
𝒜(ξ::LookupIsobar, σ::Float64) = ξ.itr(σ)

LookupIsobar(xri::Matrix{Float64}) =
    LookupIsobar(
        interpolate(
            (xri[:,1] .* mπ^2, ), # x
            xri[:,2] .+ 1im .* xri[:,3], # y
            Gridded(Linear()))) # method
# 