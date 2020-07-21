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
# 2.0 We need to define the constant of the mass we are using:
	cov = dispar.σ_we*dispar.σ2_w^0.5*dispar.σ2_e^0.5;
	d   = MvNormal([dispar.μ_w, dispar.μ_e], [dispar.σ2_w cov; cov dispar.σ2_e]);
	gg_original(θvar,evar) = pdf(d,[log(θvar),log(evar)])/(θvar*evar); #The original density function without the bounds of the distriburion
	#We now find the cummulative function to adjust the distributions:
	function integral_1(θvar) #Solving the first integral
		(val1,err1) = hcubature(x -> gg_original(θvar,x[1]), [dispar.θ_e_l],[dispar.θ_e_u]);
		return val1
	end
	(val,err)  = hcubature(x -> integral_1(x[1]), [dispar.θ_w_l],[dispar.θ_w_u]); #Solving the double integral
	cons_toMod = 1.0/val; #This is the constant to modify the distributions

# 2.1 Define the marginal distributions of the log(θ)
dist_marginal_lnθw = Normal(dispar.μ_w,dispar.σ2_w^0.5); #This is the distribution for ln(θw)
dist_marginal_lnθe = Normal(dispar.μ_e,dispar.σ2_e^0.5); #This is the distribution for ln(θe)

# 2.2 Defining the conditional distributions: for the derivative of (he) w.r.t entrepreurial ability (e)
	# 2.2.1 Distribution of θw given e/θe:
	mean_w_given_e(x_e) = dispar.μ_w + dispar.σ_we*dispar.σ2_w^0.5/dispar.σ2_e^0.5*(x_e-dispar.μ_e);
	var_w_given_e = (1.0-dispar.σ_we^2.0)*dispar.σ2_w;
	dist_w_given_e(x_e) = Normal(mean_w_given_e(x_e),var_w_given_e^0.5);
	# 2.2.2 Distribution of e/θe given θw:
	mean_e_given_w(x_θ) = dispar.μ_e + dispar.σ_we*dispar.σ2_e^0.5/dispar.σ2_w^0.5*(x_θ-dispar.μ_w);
	var_e_given_w = (1.0 - dispar.σ_we^2.0)*dispar.σ2_e;
	dist_e_given_w(x_θ) = Normal(mean_e_given_w(x_θ),var_e_given_w^0.5);

# 2.1 Define the distributions depending if there is correlation
# 2.1.1 Joint Density (gg): this is the density adjusted by the mass (is the same for both cases)
	function gg_lognorm(θvar,evar)
		# g(θ,e) = f_w|e(θ|e)*f_e(e) = f_e|w(e|θ)*f_w(θ)
		gg_original(θvar,evar)*cons_toMod;
	end #end function

	if dispar.σ_we == 0.0 #This is the case when there is no covariance between abilities
		# 2.1.2 Workers' marginal density (hw)
		function hw_lognorm_NoCov(θvar,evar)
			# h_w(θ,e) = F_e|w(e|θ)*f_w(θ)
			evar <= dispar.θ_e_l && (return 0.0)
			1.0/θvar*pdf(dist_marginal_lnθw, log(θvar))*( cdf(dist_marginal_lnθe, log(evar)) - cdf(dist_marginal_lnθe, log(dispar.θ_e_l)) )*cons_toMod;
		end #end function
		hw_lognorm = hw_lognorm_NoCov

		# 2.1.3 Entrepreneurs' marginal density (he)
		function he_lognorm_NoCov(θvar,evar)
			# h_e(θ,e) = F_w|e(θ|e)*f_e(e)
			θvar <= dispar.θ_w_l && (return 0.0)
			1.0/evar*pdf(dist_marginal_lnθe, log(evar))*( cdf(dist_marginal_lnθw, log(θvar)) - cdf(dist_marginal_lnθw, log(dispar.θ_w_l)) )*cons_toMod;
		end #end function
		he_lognorm = he_lognorm_NoCov

	else #There is a covariance in abilities
		# 2.1.2 Workers' marginal density (hw)
		function hw_lognorm_Cov(θvar,evar)
			(val,err) = hcubature(x -> gg_lognorm(θvar,x[1]), [dispar.θ_e_l],[evar]);
			return val;
			evar <= dispar.θ_e_l && (return 0.0)
		end #end function
		hw_lognorm = hw_lognorm_Cov

		# 2.1.3 Entrepreneurs' marginal density (he)
		function he_lognorm_Cov(θvar,evar)
			(val,err) = hcubature(x -> gg_lognorm(x[1],evar), [dispar.θ_w_l],[θvar]);
			return val;
			θvar <= dispar.θ_w_l && (return 0.0)
		end #end function
		he_lognorm = he_lognorm_Cov
	end #end if

# 2.1.4 Partial derivative of (he) w.r.t entrepreurial ability (e)
	function partial_he_e_lognorm(θvar,evar)
		( - 1.0/evar*(1.0 + 1.0/(1.0-dispar.σ_we^2)*(log(evar) - dispar.μ_e)/dispar.σ2_e )*he_lognorm(θvar,evar) +
		dispar.σ_we/(evar*dispar.σ2_e^0.5*dispar.σ2_w^0.5*(1.0-dispar.σ_we^2) )*(pdf(dist_marginal_lnθe, log(evar) )*( mean_w_given_e(log(evar)) -
		var_w_given_e*pdf(dist_w_given_e(log(evar)), ((log(θvar) - mean_w_given_e(log(evar)))/var_w_given_e))/cdf(dist_w_given_e(log(evar)), ( (log(θvar) - mean_w_given_e(log(evar)))/var_w_given_e) ) ) - dispar.μ_w*he(θvar,evar) ) );
	end #end function

# 3 Create distributions
	if dispar.uniform # Uniform case
		gg=gg_unif
		hw=hw_unif
		he=he_unif
		partial_he_e=partial_he_e_unif
	else # Log-Normal case
		gg = gg_lognorm	# g(θ,e) = f_w|e(θ|e)*f_e(e) = f_e|w(e|θ)*f_w(θ)
		hw = hw_lognorm # h_w(θ,e) = F_e|w(e|θ)*f_w(θ)
		he = he_lognorm # h_e(θ,e) = F_w|e(θ|e)*f_e(e)
		partial_he_e= partial_he_e_lognorm # ∂h_e/∂e
	end
	return gg, hw, he, partial_he_e
end
