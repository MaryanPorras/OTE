function marginaltaxes!(MarginalTaxes::AbstractArray{Float64,2},model::OTEmodel)
# 0 Unpack:
	#0.1 Productivities and States:
	θw  = model.states[1,:]
	e 	= model.states[2,:]
	λ	= model.prices[1]
	ω	= model.prices[2]
	#0.2 Controls:
	n	= model.controls[1,:]
	z	= model.controls[2,:]
	l	= model.controls[3,:]
	#0.3 Parameters:
	α = model.ecopar.α
	β = model.ecopar.β
	σ = model.ecopar.σ
	χ = model.ecopar.χ
	ψ = model.ecopar.ψ
	#0.4 Size of the problem:
	globalsize = model.compar.globalsize
	fullsize   = model.compar.entrepsize+globalsize-1

# 1. Define Marginal Taxes:
	#1.1 Marginal taxes:
	MarginalTaxes[1,:] = ( λ.*α.*e.*n.^(α.-1.0) .- ω )./ω #Tn'
	MarginalTaxes[2,:] = β.*z.^σ #Tc'
	MarginalTaxes[3,:] = ( θw*ω-λ*χ*l.^ψ )./(θw*ω) #Tl'

	return nothing
end

function marginaltaxes(model::OTEmodel)
	marginaltaxes = Array{Float64,2}(undef,3,model.compar.globalsize+model.compar.entrepsize-1)
	fill!(marginaltaxes,NaN)
	marginaltaxes!(marginaltaxes, model)
	return marginaltaxes
end

function taxesliabilities!(TaxesLiabilities::AbstractArray{Float64,2},model::OTEmodel,lenght_sol::Int64)
# 0 Unpack:
	# 0.1 Productivities and States:
	θw  = model.states[1,:]
	e 	= model.states[2,:]
	u 	= model.states[4,:]
	λ	= model.prices[1]
	ω	= model.prices[2]
	# 0.2 Controls:
	n	= model.controls[1,:]
	z	= model.controls[2,:]
	l	= model.controls[3,:]
	p	= model.controls[4,:]
	# 0.3 Parameters:
	α = model.ecopar.α
	β = model.ecopar.β
	σ = model.ecopar.σ
	χ = model.ecopar.χ
	ψ = model.ecopar.ψ
	# 0.4 Size of the problem:
	globalsize = model.compar.globalsize
	fullsize   = model.compar.entrepsize+globalsize-1
	solution_size = globalsize-lenght_sol+1 #The order of the last solution found in the algorithm
	# 0.5 Call the marginal taxes:
	margtaxes = marginaltaxes(model)

# 1. Define Taxes Liabilities and tax bases:
	#1.1 Tax bases:
	#Recall that for Tc we need to have the value of Tn
	TaxesLiabilities[1,:] = ω.*( n.-model.ecopar.ς ) #Tn
	TaxesLiabilities[3,:] = θw.*ω.*l #Tl

	#1.2 Taxes liabilities:
	#We need to give the initial value for the integral:
	TaxesLiabilities[4,solution_size] = 0.0 #Tn
	TaxesLiabilities[5,solution_size] = e[solution_size]*n[solution_size]^model.ecopar.α - ω/λ*(n[solution_size] - model.ecopar.ς) - model.ecopar.β/(1.0+model.ecopar.σ)*z[solution_size]^(1.0-model.ecopar.σ) - u[solution_size] #Tc
	TaxesLiabilities[6,solution_size] = θw[solution_size]*ω/λ*l[solution_size] - model.ecopar.χ/(1.0+model.ecopar.ψ)*l[solution_size]^(1.0+model.ecopar.ψ) - u[solution_size] #Tl

    #We have to get the values of the integrals:
	for i in [1,3]
		sol_int       = fill(NaN,fullsize-solution_size+1)
		base          = TaxesLiabilities[i,solution_size:fullsize]; #Tax base for the integral factor
		i == 1 ? (first = TaxesLiabilities[4,solution_size]) : (first = TaxesLiabilities[6,solution_size]) #The first element for the integral
		to_integrate  = margtaxes[i,solution_size:fullsize];
		my_integral_lb!(sol_int,to_integrate,base,first);
		i == 1 ? (TaxesLiabilities[4,solution_size:fullsize] = sol_int[:]) : (TaxesLiabilities[6,solution_size:fullsize] = sol_int[:])
	end

	# 1.3 Solve for Tc:
	# Tax base:
	TaxesLiabilities[2,:] = e.*n.^model.ecopar.α - ω.*( n.-model.ecopar.ς ) - TaxesLiabilities[4,:] - z #Tc

    #We have to get the values of the integrals:
	sol_int       = fill(NaN,fullsize-solution_size+1)
	base          = TaxesLiabilities[2,solution_size:fullsize]; #Tax base for the integral factor
	first 		  = TaxesLiabilities[5,solution_size] #The first element for the integral
	to_integrate  = margtaxes[2,solution_size:fullsize];
	my_integral_lb!(sol_int,to_integrate,base,first);
	TaxesLiabilities[5,solution_size:fullsize] = sol_int[:]

	return nothing
end

function taxesliabilities(model::OTEmodel,lenght_sol::Int64)
	taxliabilities = Array{Float64,2}(undef,6,model.compar.globalsize+model.compar.entrepsize-1)
	fill!(taxliabilities,NaN)
	taxesliabilities!(taxliabilities, model, lenght_sol)
	return taxliabilities
end
