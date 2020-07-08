# Define the static solution function
function Aiyagari(parameters,settings)

	# Unpack the parameters and settings
	@unpack (γ, α, δ, ρ, θ, σ) = parameters()
	@unpack (amin, amax, an, zmin, zmax, zn, Δ, K, Viteration, Kiteration, Vtolerance, Ktolerance, relax) = settings()

	# Construct the discretized grids
	avector = collect(range(amin, stop=amax, length=an))
	zvector = collect(range(zmin, stop=zmax, length=zn))

	# Construct the discretizes meshes
	amatrix = repeat(avector, outer=(1,zn))
	zmatrix = repeat(transpose(zvector), outer=(an,1))

	# Define the step sizes of wealth and labor efficiency
	Δa = (amax-amin)/(an-1)
	Δz = (zmax-zmin)/(zn-1)

	# Define the elements of the C matrix
	μ = (-θ*log.(zvector) .+ σ^2/2).*zvector
	σ2 = σ^2*zvector.^2
	center =  min.(μ,0)/Δz .- max.(μ,0)/Δz .- σ2/abs2(Δz)
	lower = σ2/(2*abs2(Δz)) .- min.(μ,0)/Δz
	upper = σ2/(2*abs2(Δz)) .+ max.(μ,0)/Δz

	# Construct the upper diagonal of the C matrix
	upperdiagonal = repeat([upper[1]], outer=an)
	for index in 2:zn-1
		upperdiagonal = [upperdiagonal; repeat([upper[index]], outer=an)]
	end

	# Construct the center diagonal of the C matrix
	centerdiagonal = repeat([center[1]+lower[1]], outer=an)
	for index in 2:zn-1
		centerdiagonal = [centerdiagonal; repeat([center[index]], outer=an)]
	end
	centerdiagonal = [centerdiagonal; repeat([center[zn]+upper[zn]], outer=an)]

	# Construct the lower diagonal of the C matrix
	lowerdiagonal = repeat([lower[2]], outer=an)
	for index in 3:zn
		lowerdiagonal = [lowerdiagonal; repeat([lower[index]], outer=an)]
	end

	# Construct the C matrix in the right form
	C = spdiagm(0 => centerdiagonal,
				an => upperdiagonal,
				-an => lowerdiagonal)

	# Compute prices
	r = α*K^(α-1)-δ
	w = (1-α)*K^α

	# Define the initial guess for the value function
	V = (w*zmatrix+r*amatrix).^(1-γ)/((1-γ)*ρ)

	# Define the storage arrays
	∂VF = zeros(an,zn)
	∂VB = zeros(an,zn)
	A = spdiagm(0 => ones(an*zn),
				1 => ones(an*zn-1),
				-1 => ones(an*zn-1),
				an => ones(an*(zn-1)),
				-an => ones(an*(zn-1)))

	# Define the eigenvalue correction vector
	b = fill(zero(Float64), an*zn)
	indexfix = 1
	b[indexfix] = 0.1

    # Iterate over aggregate physical capital
	for i in 1:Kiteration
		# Iterate over the value function
    	for j in 1:Viteration
		    # Forward difference
		    ∂VF[1:an-1,:] = (V[2:an,:] .- V[1:an-1,:])/Δa
		    ∂VF[an,:] = (w*zvector .+ r*amax).^(-γ) # State boundary constraint

		    # Backward difference
		    ∂VB[2:an,:] = (V[2:an,:] .- V[1:an-1,:])/Δa
		    ∂VB[1,:] = (w*zvector .+ r*amin).^(-γ) # State boundary constraint

		    # Steady state value function
		    c0 = w*zmatrix .+ r*amatrix
		    ∂V0 = c0.^(-γ)

		    # Forward drift
		    cF = ∂VF.^(-1/γ)
		    sF = w*zmatrix .+ r*amatrix .- cF
		    IF = sF .> 0

		    # Backward drift
		    cB = ∂VB.^(-1/γ)
		    sB = w*zmatrix .+ r*amatrix .- cB
		    IB = sB .< 0

		    # Flow utility
		    Vupwind = ∂VF.*IF .+ ∂VB.*IB .+ (1 .- IF .- IB).*∂V0
		    c = Vupwind.^(-1/γ)
		    u = c.^(1-γ)/(1-γ)

		    # Define the elements of the D matrix
		    center =  min.(sB,0)/Δa .- max.(sF,0)/Δa
		    lower = -min.(sB,0)/Δa
		    upper = max.(sF,0)/Δa

		    # Construct the upper diagonal of the D matrix
		    upperdiagonal = [upper[1:an-1,1]; 0]
		    for index in 2:zn-1
		        upperdiagonal = [upperdiagonal; upper[1:an-1,index]; 0]
		    end
			upperdiagonal = [upperdiagonal; upper[1:an-1,zn]]

		    # Construct the center diagonal of the D matrix
		    centerdiagonal = reshape(center, an*zn)

		    # Construct the lower diagonal of the D matrix
		    lowerdiagonal = lower[2:an,1]
		    for index in 2:zn
		        lowerdiagonal = [lowerdiagonal; 0; lower[2:an,index]]
		    end

		    # Construct the D matrix in the right form
		    D = spdiagm(0 => centerdiagonal,
		                1 => upperdiagonal,
		                -1 => lowerdiagonal)

		    # Construct the A matrix
		    A = C .+ D

		    # Construct the B matrix
		    B = (1/Δ+ρ)*sparse(I,an*zn,an*zn) .- A

		    # Reshape the value and flow utility functions
		    Vstack = reshape(V, an*zn)
		    ustack = reshape(u, an*zn)

		    # Solve the system of linear equations
		    stack = ustack .+ Vstack/Δ
		    Vstack = B\stack
		    Vnew = reshape(Vstack, (an,zn))

    		# Check convergence
    		Vdistance = maximum(abs.(Vnew .- V))
    		if Vdistance < Vtolerance
    			break
    		else
    			# Update the value function
    			V = Vnew
    		end
    	end

		# Change first column of matrix A
		for index in nzrange(A, 1)
			A.nzval[index] = zero(Float64)
		end
		A.nzval[1] = one(Float64)

		# Solve the system of linear equations
		solution = transpose(A)\b

		# Rescale the distribution function
	    correction = transpose(solution)*ones(an*zn,1)*Δa*Δz
		g = reshape(solution./correction, (an,zn));

		# Update aggregate capital
        Knew = sum(transpose(g)*avector*Δa*Δz);
		Kdistance = abs(K-Knew)
		if Kdistance < Ktolerance
			break
		else
			# Update aggregate physical capital and prices
			K = relax*K+(1-relax)*Knew
			r = α*K^(α-1)-δ
    		w = (1-α)*K^α
		end
	end

	# Return the distribution
	return g, K
end
