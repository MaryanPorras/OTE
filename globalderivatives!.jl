function globalderivatives!(du,u,pa,θ)
# 0 Unpack
	e 	= u[1]
	ϕ_e	= u[2]
	uw	= u[3]
	μ	= u[4]
	λ	= pa.prices[1]
	ω	= pa.prices[2]

# 1. Obtain densities, consider out of domain cases
	if e<pa.dispar.θ_e_l || e>pa.dispar.θ_e_u
		h_w=zero(e)
		h_e=zero(e)
		∂hw∂e=zero(e) #∂hw∂e es gg()
		∂he∂e=zero(e)	
	else
		h_e=pa.modist.he(θ, e)
		h_w=pa.modist.hw(θ, e)
		if e<pa.dispar.θ_e_l+1e-9 
			h_w=zero(e) 
			∂hw∂e=zero(e) #∂hw∂e es gg()
			∂he∂e=zero(e) 
		else
			h_w = max(h_w, 0.0)
			∂hw∂e = pa.modist.gg(θ,e) #f_θw(θ) f_θe(e)
			∂he∂e = pa.modist.partial_he_e(θ,e)	
		end
	end

# 2. Construct ss object
	ss = [θ, e, ϕ_e, uw, μ, h_e, h_w]
	pa.compar.debugbool && println("ss = ", ForwardDiff.value.(ss))

# 3. Deal with negative u
	if uw<0.0 
		du[1:6].=0.0
		du[2]=1.0
		du[4]=-1e16
		return nothing
	end

# 4. Find optimal controls
	(n, z, l, p) = globalcontrols(ss, pa.prices, pa.ecopar, pa.compar.debugbool)
	any(isnan,(n, z, l, p)) && error("Function find_derivatives gets NaN controls")

# 5. Calculate interim terms		
	h_tot= h_w + p*h_e
	Vw = pa.ecopar.utilit*uw^pa.ecopar.ϕ + θ*l*ω - λ*uw - λ*pa.ecopar.χ*l^(1.0+pa.ecopar.ψ)/(1.0+pa.ecopar.ψ);
	Ve = pa.ecopar.utilit*uw^pa.ecopar.ϕ + λ*e*n^pa.ecopar.α - λ*pa.ecopar.β*z^(1.0+pa.ecopar.σ)/(1.0+pa.ecopar.σ) - ω*(n-pa.ecopar.ς) - λ*uw;
	pa.compar.debugbool && println("(Ve*he+phi_e)*p= ", ForwardDiff.value((Ve*h_e+ϕ_e)*p))

#6. Create derivatives and return
	du[1] = p 	# ̇e' 
	du[2] = -( Vw*∂hw∂e + Ve*p*∂he∂e + λ*n^pa.ecopar.α*p*h_e) # ϕ_e'
	du[3] = pa.ecopar.χ*l^(1.0+pa.ecopar.ψ)/θ # u'
	du[4] = (λ - pa.ecopar.utilit*pa.ecopar.ϕ*uw^(pa.ecopar.ϕ-1.0))*h_tot # μ'
	du[5] = θ*l*h_w - (n-pa.ecopar.ς)*p*h_e # L'
	du[6] = ( e*n^pa.ecopar.α - pa.ecopar.β*z^(1.0+pa.ecopar.σ)/(1.0+pa.ecopar.σ) )*p*h_e - uw*h_tot - pa.ecopar.χ*l^(1.0+pa.ecopar.ψ)/(1.0+pa.ecopar.ψ)*h_w # Y'
	return nothing
end
