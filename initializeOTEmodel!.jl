function initializeOTEmodel!(model::OTEmodel, priceguess::Array{Float64,1}, finalstateguess::Array{Float64,1})
	# 0. Extract
	globalsize=model.compar.globalsize
	entrepsize=model.compar.entrepsize
	# 1. Fill θ_w values
	model.states[1,1:globalsize] = exp.(LinRange(log(model.dispar.θ_w_l+model.compar.abstol), log(model.dispar.θ_w_u), globalsize ))
	fill!(model.states[1,globalsize:end], model.dispar.θ_w_u)
	# 2. Fill e values for entrepreneurs problem
	model.states[2,globalsize:end] = exp.(LinRange(log(finalstateguess[1]), log(model.dispar.θ_e_u), entrepsize))
	# 3. Fill prices and final state 
	model.prices[:]=priceguess
	model.states[3:7, end] = finalstateguess[2:6]
end