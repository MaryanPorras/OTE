################################################################################
####																		####
####					MAIN SCRIPT FOR OPTIMAL TAXATION					####
####																		####
################################################################################

# To run inside a module, run OptTaxModule

# 0 Output switches
fig_graphs = true # Indicator to print figures (true is when we print figures).
fig_graphs && (fig_main = true) # Indicator for main figures.
fig_graphs && (fig_GloAndEnt = false) # Indicator for figures of both problems.
fig_graphs && (fig_debug = true) # Indicator for debug figures.
fig_graphs && (fig_others = false) # Indicator for debug figures.

# 1 Load packages definitions and functions
include("OTEutils.jl")

# 2 Set parameters and initialize
include("testparameters.jl")
	# This defines:
		# Parameter structures:	ecopar compar dispar
		# Initial guesses: 		priceguess finalguess
		# Specific parameters:	globalsize

# 3 Create and initialize model
model=OTEmodel(ecopar, compar, dispar)
initializeOTEmodel!(model, priceguess, finalstateguess)

# 4 Solve model
# So far this runs the inner routine.
# We need a boundary problem shooting algorithm.
lenght_sol = fullrungekutta!(model)

# 5 Generate results (DataFrames types). TODO
#		Function must receile a model and generate matrices
#Results_DF = DataFrame( hcat(transpose(model.states),transpose(model.controls)),
#						[:θw, :θe, :ϕe, :u, :μ, :L, :Y, :n, :z, :l, :p] )
#Results_DF[!, :asdf] = Results_DF[:, :n] - Results_DF[:, :z]
Results_DF = modeldata(model,lenght_sol)

# 6 Plot results. TODO
# 		Separate output vs diagnostics figures. save in Graphs Folder
fig_graphs && (fig_main && graphs_MainGlobal(Results_DF,".\\Graphs"))
fig_graphs && (fig_GloAndEnt && graphs_MainGloAndEnt(Results_DF,".\\Graphs"))
fig_graphs && (fig_debug && graphs_debug(Results_DF,".\\Graphs"))
fig_graphs && (fig_others && graphs_OtherResults(Results_DF,".\\Graphs"))
