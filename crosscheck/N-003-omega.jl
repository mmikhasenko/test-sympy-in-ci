### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 50ec2550-7256-11ec-1908-fbece8adaf13
begin
	import Pkg
	Pkg.activate(@__DIR__)
	using ThreeBodyDecay
	using Plots
	import Plots.PlotMeasures: mm
	# 
	using PlutoUI
end

# ╔═╡ 602a1012-cf64-425c-a347-c582253255dc
md"""
# Comparison of the Dalitz-plot models
"""

# ╔═╡ c5c79442-02d4-4008-a76d-b7e255cce202
theme(:wong2, grid=false, frame=:box,
	xlims=(:auto,:auto), ylim=(:auto,:auto), lab="")

# ╔═╡ 72076d51-bb4e-4eb8-91fe-66f2515fa78c
begin
	const mπ = 0.14
	const m0 = 1.04 #0.782
	const mρ = 0.77
	const Γρ = 0.15
end;

# ╔═╡ 08b66d3d-7ce6-4bf9-9971-ec799d8e5a2e
md"""
## The case  $\omega \to 3\pi$: small dalitz, tails of $\rho$
"""

# ╔═╡ f6bc93cc-6660-4021-b6d6-d75e02ad439c
const ms = ThreeBodyMasses(mπ,mπ,mπ; m0)

# ╔═╡ d23809a4-1174-4150-adb6-875ca20c567a
"""
	ϕ(σs)

Kibble function of Mandelstam invatiant variables
"""
function ϕ(σs)
	ms²=ms^2
	return	λ(λ(ms²[4],σs[1],ms²[1]),
			  λ(ms²[4],σs[2],ms²[2]),
			  λ(ms²[4],σs[3],ms²[3]))
end

# ╔═╡ 37d00c66-5784-42fc-8d05-2fbb6ee00ceb
"""
	ρ(k,σs)

Three-body phase-space factor projected to one dimension
"""
function ρ(k,σs)
	ms²=ms^2
	i,j,_ = ijk(k)
	return	sqrt(λ(ms²[4],σs[k],ms²[k])*λ(σs[k],ms²[i],ms²[j]))
end

# ╔═╡ e4147784-1a90-40c2-868a-aa41412a760b
begin
	abstract type Model end
	struct BW <: Model
		m::Float64
		Γ::Float64
	end
	amplitude(ξ::BW, σ) = ξ.m*ξ.Γ/(ξ.m^2-σ-1im*ξ.m*ξ.Γ)
end

# ╔═╡ a1c461c6-108e-4ab3-8e2b-ab5286dc26f5
function A(ξ::Model, σs)
	ϕ0=ϕ(σs)
	return ϕ0>0 ? 0.0 : sqrt(-ϕ0) * 
		  (amplitude(ξ,σs[1])+amplitude(ξ,σs[2])+amplitude(ξ,σs[3]))	
end

# ╔═╡ 336e370e-66f8-4ed2-b3b1-e67c922f5c25
md"""
## Two model: regular $\rho$ and modified $\rho$
"""

# ╔═╡ d006f277-1673-4817-a934-a1aa9a0f7b91
begin
	ξ₁ = BW(mρ,Γρ) # nominal ρ
	ξ₂ = BW(mρ, Γρ+100e-3) # modified ρ
end;

# ╔═╡ 05477d32-20cb-450b-b891-85aefeb6ca42
begin
	plot(size=(700,300), layout=grid(1,2),
		leg=:left, ylab=["|A|²" "|A|²ρ"], xlab="m(ππ) [GeV]",
		bottom_margin=5mm, left_margin=5mm)
	# 
	plot!(sp=1, m->abs2(amplitude(ξ₁,m^2)), sqrt.(lims1(ms))..., lab="model1")
	plot!(sp=1, m->abs2(amplitude(ξ₂,m^2)), sqrt.(lims1(ms))..., lab="model2")
	# 
	plot!(sp=2, m->abs2(amplitude(ξ₁,m^2))*ρ(1,m^2), sqrt.(lims1(ms))...)
	plot!(sp=2, m->abs2(amplitude(ξ₂,m^2))*ρ(1,m^2), sqrt.(lims1(ms))...)
end

# ╔═╡ 002b72b4-8276-4054-83ac-a27cf3d66eba
I(ξ::Model, σs) = abs2(A(ξ,σs))

# ╔═╡ 8b0c8a14-1eb8-4843-b0ea-a5e7eb833b37
plot(size=(700,300), layout=grid(1,2),
	plot(Base.Fix1(I, ξ₁), ms),
	plot(Base.Fix1(I, ξ₂), ms))

# ╔═╡ 245b618d-0393-479b-abc6-070939ee6421
md"""
## MC toys: likelihood ratio on pseudodata
"""

# ╔═╡ 44e36584-dc7c-49ce-bd37-258e3e92dfa3
md"""
$\mathcal{L}_i(D) = \sum_\text{D}\log |A_i(s,t)|^2$

$\Delta\mathcal{L}(D) = \mathcal{L}_1(D) - \mathcal{L}_2(D)$
"""

# ╔═╡ 74f40b1c-b447-44cb-9d7b-753b1084a03f
binneddensity(ξ::Model) = getbinned2dDensity(
	(σ1,σ3)->I(ξ, Invariants(ms;σ1,σ3)), lims1(ms), lims3(ms),  50, 50);

# ╔═╡ 81f8ef1a-f4a4-4310-be39-b1fab4de2c55
"""
	experiment(ξ₀,ξₜ,N)

Performs MC studies of the likelihood ratio.
Note: the value is shifted due to unconstrained normalization.
"""
function experiment(ξ₀,ξₜ,N)
	gds = binneddensity(ξ₀)
	pseudodata = [rand(gds) for _ in 1:N]
	σsv = [Invariants(ms; σ1, σ3) for (σ1,σ3) in pseudodata]
	Δllh = sum(σs->log(I(ξ₀, σs))-log(I(ξₜ, σs)), σsv)
	return Δllh/N
end

# ╔═╡ a535fcfb-82d0-40e4-9219-1ca51739aea4
begin
	Nev = 1000
	Nans = 200
end;

# ╔═╡ 97af4995-8b3f-4c20-a358-d435b87fc3a0
begin
	ensample1 = [experiment(ξ₁,ξ₂,Nev) for _ in 1:Nans]
	ensample2 = -[experiment(ξ₂,ξ₁,Nev) for _ in 1:Nans]
end;

# ╔═╡ 3374ea60-45ea-4a02-aa3f-60c4f1487325
begin 
	bins = range(extrema(vcat(ensample1,ensample2))..., length=30)
	plot(xlab="LLH₁-LLH₂")
	stephist!(ensample1; bins, lab="data(model1)")
	stephist!(ensample2; bins, lab="data(model2)")
end

# ╔═╡ 270bdb41-75f9-4e26-9008-5dc897136757
md"""
## Analytic expectation for toys
"""

# ╔═╡ ec2e5405-682f-44e5-924f-905e8477be06
function moment(ξ₀,ξₜ,n::Int)
	t(σs) = log(I(ξ₀,σs)) - log(I(ξₜ,σs))
	integrand(σs) = t(σs)^n * I(ξ₀,σs)
	cvalue = three_body_phase_space_integral(integrand, ms)
	return real(cvalue)
end

# ╔═╡ 308c6a64-fc93-4e69-ba25-9918c42046eb
function meansigma(ξ₀,ξₜ)
	N = moment(ξ₀,ξₜ,0)
	μN = moment(ξ₀,ξₜ,1)
	μ²N = moment(ξ₀,ξₜ,2)
	# 
	μ = μN/N
	σ = sqrt(μ²N/N - (μN/N)^2)
	(; μ, σ, N)
end

# ╔═╡ 6620abf3-dd33-46ec-974d-ce2d608567a5
begin
	μ₁,σ₁,n₁ = meansigma(ξ₁,ξ₂)
	μ₂,σ₂,n₂ = meansigma(ξ₂,ξ₁)
end;

# ╔═╡ 7d4c505e-b102-4d77-889b-54cb1e4eca15
# missing normalization
Δn = -log(n₁)+log(n₂)

# ╔═╡ 81537ba6-8d25-49a6-ba37-a20516cdf743
# Answer to the quetion 1
nσ_Q1 = (μ₂-Δn)/(σ₂/sqrt(Nev))

# ╔═╡ a17ade87-522e-41e0-96e9-b3dc589166a0
# Answer to the quetion 2
nσ_Q2 =(μ₁+μ₂)/sqrt((σ₁^2+σ₂^2)/Nev)

# ╔═╡ 7e2dc876-ac38-4144-b6be-126de67567ce
md"""
### Two questions:

__Q1__: if the nature (truth) is the model 2,
how well it is distinguished from the naive model 1 by data. The answer is __$(round(nσ_Q1, digits=1))σ__

__Q2__: we do not know what the nature is like, but how well our test statistics distinguish two models. The answer is __$(round(nσ_Q2, digits=1))σ__

The two answers differ by $\sqrt{2}$:

$n_{\sigma}^{(\text{Q2})} = \frac{2d}{\sigma_\text{eff}} = \frac{2d}{\sqrt{2}\sigma} = \sqrt{2} \frac{d}{\sigma} = \sqrt{2}\,n_{\sigma}^{(\text{Q1})}$
where
$\sigma_\text{eff} = \sqrt{\sigma_1^2+\sigma_2^2} = \sqrt{2}\sigma$
"""

# ╔═╡ 45a8e3d5-20fe-48d1-b433-3d298a3a689e
let
	σₑ(N) = sqrt((σ₁^2+σ₂^2)/N)
	# 
	answer_Q2(N) = (μ₁+μ₂)/σₑ(N) # = sqrt(N)*cosnt
	answer_Q1(N) = (μ₂-Δn)/(σ₂/sqrt(N)) # = sqrt(N)*cosnt
	# 
	plot(xlab="#events", ylab="separation in #σ", leg=:left)
	plot!(answer_Q2, 100:1_000_000, xscale=:log10, lab="Q2")
	plot!(answer_Q1, 100:1_000_000, xscale=:log10, lab="Q1")
end

# ╔═╡ Cell order:
# ╟─602a1012-cf64-425c-a347-c582253255dc
# ╠═50ec2550-7256-11ec-1908-fbece8adaf13
# ╠═c5c79442-02d4-4008-a76d-b7e255cce202
# ╠═72076d51-bb4e-4eb8-91fe-66f2515fa78c
# ╟─08b66d3d-7ce6-4bf9-9971-ec799d8e5a2e
# ╠═f6bc93cc-6660-4021-b6d6-d75e02ad439c
# ╠═d23809a4-1174-4150-adb6-875ca20c567a
# ╠═37d00c66-5784-42fc-8d05-2fbb6ee00ceb
# ╠═e4147784-1a90-40c2-868a-aa41412a760b
# ╠═a1c461c6-108e-4ab3-8e2b-ab5286dc26f5
# ╟─336e370e-66f8-4ed2-b3b1-e67c922f5c25
# ╠═d006f277-1673-4817-a934-a1aa9a0f7b91
# ╠═05477d32-20cb-450b-b891-85aefeb6ca42
# ╠═002b72b4-8276-4054-83ac-a27cf3d66eba
# ╠═8b0c8a14-1eb8-4843-b0ea-a5e7eb833b37
# ╟─245b618d-0393-479b-abc6-070939ee6421
# ╟─44e36584-dc7c-49ce-bd37-258e3e92dfa3
# ╠═74f40b1c-b447-44cb-9d7b-753b1084a03f
# ╠═81f8ef1a-f4a4-4310-be39-b1fab4de2c55
# ╠═a535fcfb-82d0-40e4-9219-1ca51739aea4
# ╠═97af4995-8b3f-4c20-a358-d435b87fc3a0
# ╠═3374ea60-45ea-4a02-aa3f-60c4f1487325
# ╟─270bdb41-75f9-4e26-9008-5dc897136757
# ╠═ec2e5405-682f-44e5-924f-905e8477be06
# ╠═308c6a64-fc93-4e69-ba25-9918c42046eb
# ╠═6620abf3-dd33-46ec-974d-ce2d608567a5
# ╠═7d4c505e-b102-4d77-889b-54cb1e4eca15
# ╠═81537ba6-8d25-49a6-ba37-a20516cdf743
# ╠═a17ade87-522e-41e0-96e9-b3dc589166a0
# ╟─7e2dc876-ac38-4144-b6be-126de67567ce
# ╠═45a8e3d5-20fe-48d1-b433-3d298a3a689e
