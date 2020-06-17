################################################################################
####																		####
####					MAIN SCRIPT FOR OPTIMAL TAXATION					####
####																		####
################################################################################

# To run inside a module, run OptTaxModule

# 0 Output switches
fig_graphs   = false # Indicator to print figures (true is when we print figures).

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
lenght_sol=fullrungekutta!(model)

# 5 Generate results (DataFrames types). TODO
#		Function must receile a model and generate matrices

# 6 Plot results. TODO
# 		Separate output vs diagnostics figures. save in Graphs Folder
