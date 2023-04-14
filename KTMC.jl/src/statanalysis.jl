
function moment(ℳ₀,ℳₜ, n::Int)
	t(σs) = log(ℐ(ℳ₀,σs)) - log(ℐ(ℳₜ,σs))
	integrand(σs) = t(σs)^n * ℐ(ℳ₀,σs)
	cvalue = three_body_phase_space_integral(integrand, masses(ℳ₀))
	return real(cvalue)
end

function meansigma(ℳ₀,ℳₜ)
	N   = moment(ℳ₀, ℳₜ, 0)
	μN  = moment(ℳ₀, ℳₜ, 1)
	μ²N = moment(ℳ₀, ℳₜ, 2)
	# 
	μ = μN/N
	σ = sqrt(μ²N/N - (μN/N)^2)
	(; μ, σ, N)
end

needN(μ₂::Real,σ₂::Real,Δn::Real; Nσ::Real) = (sqrt_N=Nσ/((μ₂-Δn)/σ₂); sqrt_N^2)

function needN(ℳ₀::Model,ℳₜ::Model; Nσ::Real)
    # 
	n₁ = moment(ℳ₀,ℳₜ,0)
	μ₂,σ₂,n₂ = meansigma(ℳₜ,ℳ₀)
    # 
    Δn = -log(n₁)+log(n₂)
    # 
    return needN(μ₂,σ₂,Δn; Nσ)
end