import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
# 
using DelimitedFiles
using KTMC
using DataFrames
using ThreeBodyDecay
using Parameters
# 
using Plots
import Plots.PlotMeasures: mm
using StatsPlots
# 
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright,
    xlim=(:auto, :auto), ylim=(:auto, :auto))
# 



df = DataFrame(folder=filter(x -> x[end] == '0', readdir(joinpath("data"))))
# 
str2m0(str) = Meta.parse(str) / 1e3
df.m0 .= str2m0.(df.folder)
sort!(df, :m0)
df = df[1:end-4, :]
# 
readinisobar(strm0::String) =
    LookupIsobar(readdlm(joinpath("data", strm0, "M_a1")))
# 
df.ξ .= readinisobar.(df.folder)
ξ₀ = readinisobar("10000")

df.N = ones(size(df.m0))
for i in 1:size(df, 1)
    @unpack m0, ξ = df[i, :]
    @show m0
    df[i, :N] = needN(
        OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ),
        OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ₀); Nσ=5)
end


df.noϕ = ones(size(df.m0))
for i in 1:size(df, 1)
    @unpack m0, ξ = df[i, :]
    @show m0
    df[i, :noϕ] = needN(
        PureSum(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ),
        PureSum(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ₀); Nσ=5)
end



df.N_just_size = ones(size(df.m0))
for i in 1:size(df, 1)
    ξ = df.ξ[12]
    m0 = df.m0[i]
    df[i, :N_just_size] = needN(
        OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ),
        OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ₀); Nσ=5)
end


df.N_just_KT = ones(size(df.m0))
for i in 1:size(df, 1)
    ξ = df.ξ[i]
    m0 = df.m0[20]
    df[i, :N_just_KT] = needN(
        OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ),
        OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ₀); Nσ=5)
end


let
    @df df plot(:m0, :N, yscale=:log10, m=(4, :black), lab="1++",
        xlab="mass of the system [GeV]", ylab="# events required for 5σ")
    @df df plot!(:m0, :N_just_size, yscale=:log10, m=(3, :o), mc=2, lab="just size",
        xlab="mass of the system [GeV]", ylab="# events required for 5σ")
    @df df plot!(:m0, :N_just_KT, yscale=:log10, m=(3, :o), mc=3, lab="just KT",
        xlab="mass of the system [GeV]", ylab="# events required for 5σ")
end
savefig(joinpath("plots", "1--KT.pdf"))


#                      _|                      _|                      
#  _|      _|      _|        _|_|_|    _|_|_|  _|    _|_|      _|_|_|  
#  _|      _|      _|  _|  _|    _|  _|    _|  _|  _|_|_|_|  _|_|      
#    _|  _|  _|  _|    _|  _|    _|  _|    _|  _|  _|            _|_|  
#      _|      _|      _|    _|_|_|    _|_|_|  _|    _|_|_|  _|_|_|    
#                                _|        _|                          
#                            _|_|      _|_|                            

shapeeff = DataFrame(m0=range(0.6, 4.0, length=100))
for iξ = 10:2:20
    ξ = df.ξ[iξ]
    m0_KT = df.folder[iξ]
    shapeeff[:, "N_m0_"*string(m0_KT)] .= ones(size(shapeeff.m0))
    # 
    for i in 1:size(shapeeff, 1)
        m0 = shapeeff.m0[i]
        shapeeff[i, "N_m0_"*string(m0_KT)] = needN(
            OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ),
            OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ₀); Nσ=5)
    end
end

let
    plot(xlab="mass of the system [GeV]", ylab="# events required for 5σ")
    @df shapeeff plot!(:m0, :N_m0_0900, yscale=:log10, lab="mKT=0.9 GeV")
    @df shapeeff plot!(:m0, :N_m0_1000, yscale=:log10, lab="mKT=1.0 GeV")
    @df shapeeff plot!(:m0, :N_m0_1100, yscale=:log10, lab="mKT=1.1 GeV")
    @df shapeeff plot!(:m0, :N_m0_1500, yscale=:log10, lab="mKT=1.5 GeV")
    @df shapeeff plot!(:m0, :N_m0_2000, yscale=:log10, lab="mKT=2.0 GeV")
end
savefig(joinpath("plots", "wiggles.pdf"))




let i = 12
    @unpack ξ, m0 = df[i, :]
    ℳ₀ = OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ₀)
    ℳ = OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ)
    #
    plot(layout=grid(1, 3), size=(1200, 300),
        plot(ℳ₀, colorbar=true, title="full (m₀=$(m0) GeV)"),
        plot(σs -> log(ℐ(ℳ₀, σs) / ℐ(ℳ, σs)), KTMC.masses(ℳ), colorbar=true, title="log(I₀/I)"),
        plot(σs -> ℐ(ℳ₀, σs) * log(ℐ(ℳ₀, σs) / ℐ(ℳ, σs)), KTMC.masses(ℳ), colorbar=true, title="I₀ log(I₀/I)"))
    # 
    # savefig(joinpath("plots","dalitzNLL_m0$(df.folder[i]).pdf"))
end



nointer(ℳ, σs) = -Kibble(σs, KTMC.masses(ℳ)^2) *
                 (ℐ(ℳ.ξ, σs[1]) + ℐ(ℳ.ξ, σs[2]) + ℐ(ℳ.ξ, σs[3]))
# 
interference(ℳ, σs) = ℐ(ℳ, σs) - nointer(ℳ, σs)
# 
let i = 24
    @unpack ξ, m0 = df[i, :]
    ℳ₀ = OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ₀)
    #
    p1 = plot(layout=grid(1, 3), size=(1200, 300),
        plot(ℳ₀, colorbar=true, title="full"),
        plot(Base.Fix1(nointer, ℳ₀), KTMC.masses(ℳ₀), colorbar=true, title="no interference"),
        plot(Base.Fix1(interference, ℳ₀), KTMC.masses(ℳ₀), colorbar=true, title="just interference"))
    #
    ℳ = OneMinusMinus(ThreeBodyMasses(mπ, mπ, mπ; m0), ξ)
    p2 = plot(layout=grid(1, 3), size=(1200, 300),
        plot(σs -> log(ℐ(ℳ₀, σs) / ℐ(ℳ, σs)), KTMC.masses(ℳ), colorbar=true, title="full"),
        plot(Base.Fix1(nointer, ℳ), KTMC.masses(ℳ), colorbar=true, title="no interference"),
        plot(Base.Fix1(interference, ℳ), KTMC.masses(ℳ), colorbar=true, title="just interference"))
    #
    plot(p1, p2, layout=grid(2, 1), size=(1200, 650))
end

