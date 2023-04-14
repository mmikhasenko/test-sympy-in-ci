import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using KTMC
using ThreeBodyDecay

using Plots
import Plots.PlotMeasures: mm
using LaTeXStrings
# 
theme(:wong2, grid=false, frame=:box,
	xlims=(:auto,:auto), ylim=(:auto,:auto), lab="")
# 

criticalpoints(m,mR) = Dict(
    L"m_R+m" => mR+m,
    L"\sqrt{2m_R^2+m^2}" => sqrt(2mR^2+m^2),
    L"\sqrt{3m_R^2-3m^2}" => sqrt(3mR^2-3m^2),
    L"(m_R^2-m^2) / m" => (mR^2-m^2) / m)

let m = 0.3, mR=0.7
    plot(layout=grid(2,2), size=(500,550))
    cp = criticalpoints(m,mR)
    for (sp,(t,m0)) in enumerate(cp)
        plot!(PureSum(ThreeBodyMasses(m,m,m; m0), KTMC.BW(mR,0.1Γρ)); sp, ticks=false, title=L"m_0="*t)
    end
    plot!()
end

savefig(joinpath("plots", "topologicaltransitions.pdf"))


modelvsphsp = let
    xv = range(0.8, 2.5, length=100)
    yv = [meansigma(
        PhaseSpace(ThreeBodyMasses(mπ,mπ,mπ; m0)),
        OneMinusMinus(ThreeBodyMasses(mπ,mπ,mπ; m0), KTMC.BW(mρ, 0.1Γρ))) for m0 in xv]
    (; xv, yv)
end

let
    plot(modelvsphsp.xv, getproperty.(modelvsphsp.yv, :μ),
        xlab = L"\mathrm{mass\,\,of\,\,the\,\,system\,\,[MeV]}",
        ylab=L"\int \log(I(s,t))\,\mathrm{d} s\,\mathrm{d} t ")
    # 
    cp = criticalpoints(mπ,mρ)
    for (t,m0) in cp
        vline!([m0], lab=t)
    end
    plot!(xlim=(0.8,2.5))
end
savefig(joinpath("plots", "purelogarithm.pdf"))

plot(
    [plot(ThreeBodyMasses(mπ,mπ,mπ; m0),
        σs->log(ℐ(
            OneMinusMinus(
                ThreeBodyMasses(mπ,mπ,mπ; m0),
                KTMC.BW(mρ,0.1Γρ)), σs)),
        title=latexstring("m_0=$(m0)\\,\\mathrm{GeV}"))
        for m0 in 1.3:0.1:1.8]...)
savefig(joinpath("plots", "interferencering.pdf"))
