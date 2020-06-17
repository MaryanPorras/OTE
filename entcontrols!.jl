function entcontrols!(controls::AbstractArray{T,1}, ss::Array{T,1}, prices::Array{Float64,1}, pa::EconomicParameters, verbosebool::Bool) where T<:Real
	# INPUT: control vector, state vector and parameters
	# MODIFIES: control vector
# 1.0 Preallocating
	corner_nn::T   =NaN
	corner_zz::T   =NaN
	corner_hamiltonian::T=NaN
	max_hamiltonian::T=-Inf
# 1.1 Extracting
	e::T=ss[1]
	u::T=ss[2]
	μ::T=ss[3]
	h_e::T =ss[4]
	λ::Float64=prices[1]
	ω::Float64=prices[2]

# 1.2 Initialize constants
	z_u::T	= (1.0/pa.β)^(1.0/pa.σ)*(1.0 - 1e-10) #Max possible evasion.
	ln_max::T	= 1e10
	n_full_info::T  = (λ*pa.α*e/ω)^(1.0/(1.0-pa.α))
	n_u::T		= min((-λ*z_u*h_e/μ/pa.σ)^(1.0/pa.α), n_full_info)
	pre_tax_profits_at_nmin::T = λ*e*pa.ς^pa.α # - ω*pa.ς

# 1.3 Define aux functions:
	z_opt(nvar)     = - pa.σ*μ*nvar^pa.α/(λ*h_e);
	nfoc(nvar,zvar) = λ*pa.α*e*nvar^(pa.α-1.0)*h_e - ω*h_e + pa.α*μ*nvar^(pa.α-1.0)*(1.0-pa.β*zvar^pa.σ);
	nfoc_at_zopt(nvar) = nfoc(nvar,z_opt(nvar))
	objective(nvar,zvar) = ( pa.utilit*u^pa.ϕ*h_e 
							+ λ*( e*nvar^pa.α - pa.β*zvar^(1.0+pa.σ)/(1.0+pa.σ) - u)*h_e
							 - ω*nvar*h_e + μ*nvar^pa.α*(1.0-pa.β*zvar^pa.σ) )




	n_opt(zvar)         = pa.α/ω/(1.0-pa.α)*(λ*pa.β/(1.0+pa.σ)*zvar^(1.0+pa.σ)-A_cons)
	zfoc(nvar, zvar)    = λ*e*nvar^pa.α - ω*nvar - λ/pa.σ*zvar*(1.0-pa.β/(1.0+pa.σ)*zvar^pa.σ) + A_cons
	zfoc_at_nopt(zvar)  = zfoc(n_opt(zvar), zvar)
	zfoc_at_nmin(zvar)  = zfoc(pa.ς, zvar)

# 2 Define cases and solve

# 2.0 Check μ=0.0 and solve
	if abs(μ)<1e-12
		controls[1]=n_full_info
		controls[2]=0.0
		verbosebool && println("n = ", ForwardDiff.value(controls[1]), ",  z = ", ForwardDiff.value(controls[2]))
		return nothing
	end

# 2.1 Evaluate potential interior solution:
	# Using bisection method:
	if nfoc_at_zopt(pa.ς)*nfoc_at_zopt(n_u)<0.0
		verbosebool && println("Case 1: Interior n")
		controls[1] = find_zero(nfoc_at_zopt, (pa.ς, n_u), Bisection())
		controls[2] = z_opt(controls[1])
		max_hamiltonian = objective(controls[1],controls[2])
	end

# 2.2 Evaluate low n corner solution:
	corner_nn = pa.ς
	corner_zz = z_opt(corner_nn)
	corner_hamiltonian = objective(corner_nn,corner_zz)

	# 2.2.1 Keep best solution
	if corner_hamiltonian>max_hamiltonian
		verbosebool && println("Case 2: Corner n")
			controls[1]=corner_nn
			controls[2]=corner_zz
			max_hamiltonian=corner_hamiltonian
	end

# 2.3 Case 3: n at maximum value
	corner_nn = n_u
	corner_zz = z_opt(corner_nn)
	corner_hamiltonian = objective(corner_nn,corner_zz)
	if corner_hamiltonian > max_hamiltonian
		controls[1]   = corner_nn
		controls[2]   = corner_zz
		max_hamiltonian = corner_hamiltonian
	end

# 3 Final Output:
	verbosebool && println("n = ", ForwardDiff.value(controls[1]), ",  z = ", ForwardDiff.value(controls[2]))
end

function entcontrols(ss::Array{T,1}, prices::Array{Float64,1}, pa, verbosebool::Bool=false) where T<:Real
	tempcont=Array{T}(undef,2)
	fill!(tempcont,NaN)
	entcontrols!(tempcont, ss, prices, pa, verbosebool)
	return Tuple(tempcont)
end
