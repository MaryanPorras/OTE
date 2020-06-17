function createdistributions(dispar::DistributionParameters)
# 1 Bivariate uniform abilities 
# 1.1 Joint Density (gg)
	function gg_unif(θvar,evar)	
		# g(θ,e) = f_w|e(θ|e)*f_e(e) = f_e|w(e|θ)*f_w(θ)
		1.0/(dispar.θ_w_u-dispar.θ_w_l)/(dispar.θ_e_u-dispar.θ_e_l)	
	end
# 1.2 Workers' marginal density (hw)
	function hw_unif(θvar,evar)
		# h_e(θ,e) = F_w|e(θ|e)*f_e(e)
		(evar-dispar.θ_e_l)/(dispar.θ_e_u-dispar.θ_e_l)/(dispar.θ_w_u-dispar.θ_w_l) 
	end
# 1.3 Entrepreneurs' marginal density (he)
	function he_unif(θvar,evar)
		# h_w(θ,e) = F_e|w(e|θ)*f_w(θ)
		(θvar-dispar.θ_w_l)/(dispar.θ_w_u-dispar.θ_w_l)/(dispar.θ_e_u-dispar.θ_e_l) 
	end
# 1.4 Partial derivative of (he) w.r.t entrepreurial ability (e)
	function partial_he_e_unif(θvar,evar)
		0.0
	end

# 2 Bivariate log-normal abilities

# 3 Create distributions
	if dispar.uniform # Uniform case
		gg=gg_unif
		hw=hw_unif
		he=he_unif
		partial_he_e=partial_he_e_unif
	else # Log-Normal case
		hw= NaN # h_w(θ,e) = F_e|w(e|θ)*f_w(θ)
		he= NaN # h_e(θ,e) = F_w|e(θ|e)*f_e(e)
		gg= NaN	# g(θ,e) = f_w|e(θ|e)*f_e(e) = f_e|w(e|θ)*f_w(θ)
		partial_he_e= NaN	# ∂h_e/∂e
	end
	return gg, hw, he, partial_he_e
end

