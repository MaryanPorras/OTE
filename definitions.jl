# Notes and reminders
# using Distributions
# using FastGaussQuadrature

# 1. Parameters
# 1.1 Economic parameters
struct EconomicParameters
	# 1.1.1 Workers' parameters
	χ::Float64 	# Scale parameter for labor supply
	ψ::Float64 	# Elasticity of labor supply
	# 1.1.2 Entrepreneurs' parameters
	α::Float64 	# Returns to scale/labor elasticity of output
	β::Float64	# Scale parameter for evasion
	σ::Float64 	# Elasticity of evasion
	ς::Float64 	# Entrepreneurs' self-provision of labor hours
	# 1.1.3 Labor informality parameters
	informal::Bool	# Indicator of active-inactive informality
	κ::Float64	# Workers' informal supply scale parameter
	ρ::Float64	# Workers' informal supply elasticity
	δ::Float64	# Entrepreners' informal demand scale parameter
	γ::Float64	# Entrepreneurs' informal demand elasticity
	# 1.1.4 Planner's parameters
	utilit::Bool	# Indicator of utilitarian planner
	ϕ::Float64 	# Concave utilitarian parameter
	G::Float64 	# Expenditure needs
end

#################################################################
#																#
#	NEED TO BETTER DEFINE DISTRIBUTION PARAMETER CONSTRUCTORS	#
#																#
#################################################################

# 1.2 Distribution parameters
struct DistributionParameters
	# 1.2.1 Type of distribution
	uniform::Bool	# Indicator of uniform distributions.
	# 1.2.2 Moments
	μ_w::Float64	# Mean of worker's ability
	μ_e::Float64	# Mean of entrepreneur's ability
	σ2_w::Float64	# Variance of worker's ability
	σ2_e::Float64	# Varaince of entrepreneur's ability
	σ_we::Float64	# Covariance of abilities.
	# 1.2.3 Support
	θ_w_l::Float64	# Lower bound of worker's ability
	θ_w_u::Float64	# Upper bound of worker's ability
	θ_e_l::Float64	# Lower bound of entrepreneurs's ability
	θ_e_u::Float64	# Upper bound of entrepreneurs's ability
end

# 1.3 Computation parameters
struct ComputationParameters
	globalsize::Int64
	entrepsize::Int64
	abstol::Float64
	reltol::Float64
	verbosebool::Bool
	debugbool::Bool
	alg::OrdinaryDiffEqAlgorithm
end

# 2. Distrebutions and composite paramenters

# 2.1 Model distribution functions
	# 2.1.1 Load funciton that creates distributions from parameters
include("createdistributions.jl")
	#2.1.2 Define struct and constructor
struct ModelDistributions
	gg::Function
	hw::Function
	he::Function
	partial_he_e::Function
	function ModelDistributions(dispar::DistributionParameters)
		(gg, hw, he, partial_he_e)=createdistributions(dispar)
		new(gg, hw, he, partial_he_e)
	end
end

# 2.2 # Wrapper to pass parameters to Runge-Kutta solvers
struct RKParameters 
	ecopar::EconomicParameters
	compar::ComputationParameters
	dispar::DistributionParameters
	modist::ModelDistributions
	prices::Array{Float64,1}	
end

# 3. Model definition
struct OTEmodel
	# 3.1 Elements
	ecopar::EconomicParameters
	compar::ComputationParameters
	dispar::DistributionParameters
	modist::ModelDistributions
	states::Array{Float64,2} 	
	controls::Array{Float64,2}	
	prices::Array{Float64,1}	# Constant states
	# 3.2 Constructor
	function OTEmodel(ecopar::EconomicParameters, compar::ComputationParameters, dispar::DistributionParameters)
		modist	=ModelDistributions(dispar)
		states	=Array{Float64}(undef, 7, compar.globalsize+compar.entrepsize-1)	# Order: θ, e, ϕ_e, u, μ, L, Y
		controls=Array{Float64}(undef, 4, compar.globalsize+compar.entrepsize-1)	# Order: n, z, l, p
		prices	=Array{Float64}(undef, 2)											# Order: λ, ω
		fill!(states, NaN)
		fill!(controls, NaN)
		fill!(prices, NaN)  
		new(ecopar, compar, dispar, modist, states, controls, prices)
	end	
end

