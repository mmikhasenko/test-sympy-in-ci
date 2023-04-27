import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using KTMC
using KTMC.ThreeBodyDecay

using Plots
import Plots.PlotMeasures: mm

# 
theme(:wong2, grid=false, frame=:box,
	xlims=(:auto,:auto), ylim=(:auto,:auto), lab="")
# 


const ms = ThreeBodyMasses(mπ,mπ,mπ; m0=mϕ)

ℳ₁ = OneMinusMinus(ms, KTMC.BW(mρ,Γρ))
ℳ₂ = OneMinusMinus(ms, KTMC.BW(mρ, Γρ+100e-3))

needN(ℳ₁,ℳ₂; Nσ=5)

plot(δm->needN(
    OneMinusMinus(ms, KTMC.BW(mρ,Γρ)),
    OneMinusMinus(ms, KTMC.BW(mρ+δm/1e3, Γρ)); Nσ=5), range(-100,100, length=50),
    yscale=:log10)

#
const mρ′ = mρ+0.2

m0deps = let
    xv = range(mω, 2.0, length=100)
    yv = [needN(
        OneMinusMinus(ThreeBodyMasses(mπ,mπ,mπ; m0), KTMC.BW(mρ,Γρ)),
        OneMinusMinus(ThreeBodyMasses(mπ,mπ,mπ; m0), KTMC.BW(mρ′, Γρ)); Nσ=5) for m0 in xv]
    (; xv, yv)
end
# 

let
    plot(m0deps.xv, m0deps.yv, yscale=:log10,
        xlab = "mass of the system [MeV]", ylab="# events required for 5σ")
    vspan!(mρ .+ Γρ/2 .* [-1, 1], lab="m(ρ)", lw=0.0, α=0.2)
    vspan!(mρ′ .+ Γρ/2 .* [-1, 1], lab="m'(ρ)", lw=0.0, α=0.2)
end

plot(OneMinusMinus(ThreeBodyMasses(mπ,mπ,mπ; m0=1.65), KTMC.BW(mρ,Γρ)))