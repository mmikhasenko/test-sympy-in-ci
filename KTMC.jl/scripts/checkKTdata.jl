using DelimitedFiles
using QuadGK

using KTMC

using Plots
import Plots.PlotMeasures: mm

# 
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))
# 


const path = joinpath("data","0450","M_a1")

str2m0(str) = Meta.parse(str) / 1e3
const allm0 = sort(filter(x->x[end]=='0', readdir(joinpath("data"))), by=str2m0)

readinisobar(strm0::String) =
    LookupIsobar(readdlm(joinpath("data",strm0,"M_a1")))

ξ₀ = readinisobar(allm0[end])

ξs = [
    "Omnes"  => readinisobar(allm0[end]),
    "ω"      => readinisobar("0850"),
    "ϕ"      => readinisobar("1050"),
    "1.5GeV" => readinisobar("1500"),
    "J/ψ"    => readinisobar("3000")
]


arg(z;ϕ=π/2) = atan(imag(z*cis(-ϕ)), real(z*cis(-ϕ)))+ϕ

let emax = 2.1
    plot(layout=grid(1,4), size=(1600,350), xlab="m(ππ) (GeV)", bottom_margin=8mm)
    for (lab, ξ) in ξs
        n = sqrt(quadgk(σ->ℐ(ξ, σ), (2mπ)^2, emax^2)[1])
        plot!(sp=1, e->ℛℯ(ξ, e^2)/n, 2mπ, emax; lab, title="Real part")
        plot!(sp=2, e->ℐ𝓂(ξ, e^2)/n, 2mπ, emax; lab, title="Imaginary part")
        plot!(sp=3, e->ℐ(ξ, e^2)/n^2, 2mπ, emax; lab, title="Squared amplitude")
        plot!(sp=4, e->arg(𝒜(ξ, e^2)), 2mπ, emax; lab, leg=:bottom, title="Phase")
    end
    plot!()
end