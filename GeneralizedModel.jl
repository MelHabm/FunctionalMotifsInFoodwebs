module GeneralizedModel

using Distributions,
UUIDs,
LinearAlgebra,
StatsBase

export gen_J, assign_lifespans, determine_trophic_levels

# -------------------------------------------------------------------------------------------------------------------
abstract type AbstractSpecies end

struct Species <: AbstractSpecies
    np::Float64         # niche position
    r::Float64          # feeding range
    c_fr::Float64       # center of the feeding range
    Tpk::Float64        # temperature peak
    id::UUID            # unique species ID
end

abstract type AbstractCommunity end

"""
    Community

Type containing parameters defining a community including its interaction matrix and a vector of species within.
"""
struct Community <: AbstractCommunity
    N::Int
    A::Matrix{Float64}
    sp::Vector{Species}
    ids::Vector{UUID}
    np::Vector{Float64}
    T::Float64
    R::Float64
end

struct StructuralParameters
    N::Int
    A::Matrix{Float64}
    α::Vector{Float64}
    β::Matrix{Float64}
    χ::Matrix{Float64}
    ρ::Vector{Float64}
    ρ̃::Vector{Float64}
    σ::Vector{Float64}
    σ̃::Vector{Float64}
end

struct ExponentialParameters
    γ::Vector{Float64}
    λ::Matrix{Float64}
    μ::Vector{Float64}
    ϕ::Vector{Float64}
    ψ::Vector{Float64}
end

"""
    GeneralisedParameters

Type to contain the parameters for the normalised model.
"""
struct GeneralisedParameters
    s::StructuralParameters
    e::ExponentialParameters
end

"""
gen_J()

Function that generates a Jacobian matrix of a stable system based on a Generalized Foodweb Model (based on the model in 
Gross and Feudel 2006)
"""
function gen_J(A; R::Float64 = 42.0)
    N = size(A, 1);

    # Network topology is based on a niche model. 
    # Since the example network has no niche values needed for determining 
    # the time scale, first calculate niche like values based on the trophic level of a species
    trophic_levels = determine_trophic_levels(A)  
    niche_like_values = trophic_to_niche(trophic_levels)

    # Create an empty Jacobian with the same size as the interaction matrix
    J = zeros(size(A));

    # Calculate alpha based on metabolic rate scales from the niche-like values
    α = R .^ (0.25 .* niche_like_values)

    # Fill the entries in the Jacobian matrix using a generalized food web model 
    # (for more information about the model read Gross & Feudel 2006)

    # Create stable system 
    # A system is stable if the real part of the largest eigenvalue is negative
    while maximum(real(eigvals(J))) >= 0 

        # pick scale and exponent parameters
        ϕ = rand(Uniform(0.0, 1.0), N)
        ρ = 1 .* (sum(A, dims = 2) .!= 0)[:]   
        ρ_tilde = 1 .- ρ
        γ = rand(Uniform(0.5, 1.5), N)
        λ = ones(N, N)  
        ψ = rand(Uniform(0.5,1.5), N)
        σ = 1 .* (sum(A, dims = 1) .!= 0)[:]
        σ_tilde = 1 .- σ
        μ = rand(Uniform(1.0, 2.0), N)

        β = A ./ norm.(eachcol(A), 1)'  # creates NaN when divisor is 0
        χ = A ./ norm.(eachrow(A), 1)

        β[isnan.(β)] .= 0.0;
        χ[isnan.(χ)] .= 0.0;

        for i = 1:N # loop over rows
            for j = 1:N # loop over columns

                if i == j # intraspecific dynamic
                    J[i,j] = ρ_tilde[i] * ϕ[i] +                                     # primary production
                             ρ[i] * (γ[i] * χ[i,i] * λ[i,i] + ψ[i]) -                # predation gain
                             σ_tilde[i] * μ[i]                                       # natural mortality
                    for k = 1:N
                        J[i,j] -= σ[i] * β[k,i] * λ[k,i] * ((γ[k] - 1) * χ[k,i] + 1) # predation loss
                    end
                else      # interspecific dynamic
                    J[i,j] = 0

                    if χ[i,j] != 0
                        J[i,j] = ρ[i] * γ[i] * χ[i,j] * λ[i,j] 
                    end

                    if β[j,i] != 0
                        J[i,j] -= σ[i] * β[j,i] * ψ[j]
                    end
    
                    for k = 1:N
                        if (β[k,i] != 0) && (χ[k,j] != 0)
                            J[i,j] -= σ[i] * (β[k,i] * λ[k,j] * (γ[k] - 1) * χ[k,j])
                        end
                    end
                end
    
                J[i, j] *= α[i]
            end
        end
    end

    return J
end

# Assign lifespans based on trophic levels
function assign_lifespans(trophic_levels)
    base_lifespan = 1.0       # Lifespan for the basal trophic level
    lifespan_increment = 3.0  # Increment for each higher trophic level

    lifespans = base_lifespan .+ (trophic_levels .- 1) .* lifespan_increment
    return lifespans
end

# Determine the trophic level
function determine_trophic_levels(Interactions)
    # Number of species
    n = size(Interactions, 1)

    # Initialize trophic levels with ones
    trophic_levels = ones(Float64, n)

    # Iteratively update trophic levels until convergence
    max_iterations = 1000
    tolerance = 1e-6

    for _ in 1:max_iterations
        new_trophic_levels = ones(Float64, n)
        for i in 1:n
            prey_indices = findall(x -> x == 1, Interactions[i, :])
            if !isempty(prey_indices)
                new_trophic_levels[i] += mean(trophic_levels[prey_indices])
            end
        end
        
        # Check for convergence
        if norm(new_trophic_levels - trophic_levels) < tolerance
            break
        end
        
        trophic_levels = new_trophic_levels
    end
    
    return trophic_levels
end

"""

    structural_parameters(N, A, np, R)

Define structural parameters including

Scale parameters = rate of biomass flow in the steady state
α   : characteristic times scales
Relative contributions:
ρ   : fraction of growth rate by predation
ρ̃   : fraction of growth by primary production
σ   : fraction of loss because of predation
σ̃   : fraction of loss because of mortality
β   : per capita loss rate because of predation by another species
χ   : relative contribution to total prey of another species
"""

function structural_parameters(N, A, np, R)
    α = (R .^ 0.25np)

    β = A ./ norm.(eachcol(A), 1)'
    χ = A ./ norm.(eachrow(A), 1)

    β[isnan.(β)] .= 0.0
    χ[isnan.(χ)] .= 0.0

    ρ = 1 .* (sum(A, dims = 2) .!= 0)[:]
    ρ̃  = 1 .- ρ
    
    σ = 1 .* (sum(A, dims = 1) .!= 0)[:]
    σ̃  = 1 .- σ

    return StructuralParameters(N, A, α, β, χ, ρ, ρ̃, σ, σ̃)
end

"""

    random_parameters(N::Int64, M::Int64)

Create a set of exponent parameters for the Generalized Foodweb Model.
They describe the nonlinearity of the ecological processes in the steady state.
"""

function random_parameters(N::Int64, M::Int64)
    γ = rand(Uniform(0.5, 1.5), N, M)       # non-linearity of predation rate with respect to prey Density
    λ = ones(N, N)                          # non-linearity of the contribution to the diet of another species
    μ = rand(Uniform(1.0, 2.0), N, M)       # non-linearity of mortality
    ϕ = rand(Uniform(0.0, 1.0), N, M)       # non-linearity of the primary production
    ψ = rand(Uniform(0.5,1.5), N, M)        # non-linearity of the response function's predator dependence

    return [ExponentialParameters(γ[:,i], λ, μ[:,i], ϕ[:,i], ψ[:,i]) for i = 1:M]
end

"""
    generalised_jacobian!(J, s::StructuralParameters, e::ExponentialParameters)

Constructs the Jaccobian matrix following a Generalized Foodweb Model. 
Input:
N : community size
s : struct with structural parameters
e : struct with exponent parameters
"""
function generalised_jacobian!(J, s::StructuralParameters, e::ExponentialParameters)

    # Generate the Jacobian
    for i = 1:s.N
        for j = 1:s.N
            if i == j   # intraspecific dynamics
                J[i, j] = s.ρ̃[i] * e.ϕ[i] +                                 # Production 
                         s.ρ[i] * (e.γ[i] * s.χ[i,i] * e.λ[i,i] + e.ψ[i]) - # Growth through predation 
                         s.σ̃[i] * e.μ[i]                                    # Mortality
                for k = 1:s.N
                    # Substract predation losses
                    J[i, j] -= s.σ[i] * s.β[k, i] * e.λ[k, i] * ((e.γ[k] - 1) * s.χ[k, i] + 1) 
                end
            else        # interspecific dynamics
                J[i,j] = 0

                if s.χ[i, j] != 0
                    J[i, j] = s.ρ[i] * e.γ[i] * s.χ[i,j] * e.λ[i,j] 
                end

                if s.β[j, i] != 0
                    J[i, j] -= s.σ[i] * s.β[j,i] * e.ψ[j]
                end

                for k = 1:s.N
                    if (s.β[k,i] != 0) && (s.χ[k,j] != 0)
                        J[i,j] -= s.σ[i] * (s.β[k,i] * e.λ[k,j] * (e.γ[k] - 1) * s.χ[k,j])
                    end
                end
            end

            J[i,j] *= s.α[i]
        end
    end
end

"""
    generalised_jacobian(N::Int, s::StructuralParameters, e::ExponentialParameters)

Generate the Jacobian matrix.

Input:
N : community size
s : struct with structural parameters
e : struct with exponent parameters
"""
function generalised_jacobian(N::Int, s::StructuralParameters, e::ExponentialParameters)
    J = zeros(N, N)
    generalised_jacobian!(J, s, e)
    return J 
end

"""

    trophic_to_niche(trophic_levels)

Function to get an approximation of niche values based on a given network structure.
""" 
function trophic_to_niche(trophic_levels)
    # Normalize trophic levels to range [0, 1]
    min_level = minimum(trophic_levels)
    max_level = maximum(trophic_levels)
    
    # Avoid division by zero if all trophic levels are equal
    if max_level == min_level
        return zeros(length(trophic_levels))
    else
        niche_values = (trophic_levels .- min_level) ./ (max_level - min_level)

        # Add noise to come closer to the niche values used in the niche model
        niche_values = clamp.(niche_values .+ randn(length(niche_values)) * 0.05, 0, 1)
        return niche_values
    end
end

end # end module