

# 1. Load packages

# 1.1 Distributions and non-linear solvers
using Statistics
using Roots

# 1.2 Differential equations 
using DifferentialEquations
using ForwardDiff
# 1.2.1 Extend prevfloat for duals
	Base.nextfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(nextfloat(d.value), d.partials) 
	Base.prevfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(prevfloat(d.value), d.partials)

# 1.3 Saving output
using DataFrames
using CSV
# 1.4 Make figures
fig_graphs && using PyPlot

# 2. Type definitions
include("definitions.jl")

# 3. Solving Functions
# 3.1 Initialization
include("initializeOTEmodel!.jl")
# 3.2 Runge-Kutta solver functions
include("fullrungekutta!.jl")
include("entderivatives!.jl")
include("entcontrols!.jl")
include("globalderivatives!.jl")
include("globalcontrols!.jl")


# 4 Results functions

	#The file for the function that computes everything:
	# include("ProblemFunction.jl")
	#Marginal taxes:
	# include("marginal_taxes.jl")
	# Propositions:
	# include("Propositions.jl")
	# include("Integrals.jl")

