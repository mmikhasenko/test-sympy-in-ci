using Test
using KTMC
using KTMC.ThreeBodyDecay



ms = ThreeBodyMasses(mπ, mπ, mπ; m0=mϕ)
σs = randomPoint(ms)

ℳ₁ = OneMinusMinus(ms, KTMC.BW(mρ, Γρ))
ℳ₂ = OneMinusMinus(ms, KTMC.BW(mρ, Γρ + 100e-3))
Φ = PhaseSpace(ms)

@testset "Model tests" begin
    @test 𝒜(ℳ₁, σs) isa Complex
    @test 𝒜(ℳ₂, σs) isa Complex
    @test 𝒜(Φ, σs) isa Complex
    # 
    @test 𝒜(Φ, σs) == 1.0
end

@testset "need N" begin
    @test needN(ℳ₁, ℳ₂; Nσ=5) ≈ 5138.914212258855
end
