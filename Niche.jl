module Niche

using Distributions,
UUIDs,
LinearAlgebra,
Graphs

"""
Foodweb generator based on the ecological niche model.
The following functions are used to:
- create species with certain traits (niche parameters)
- build a community as a network of species with predator-prey links
- add/remove/move species between communities
- compare communities for equality
- avoid cases like cannibalism or mutual predation

Code is based on code provided by Tom Clegg (https://github.com/CleggTom/fw_reactivity)
"""

"""
Create custom data types as struct

Species

Type to store niche model parameters for a single species, its thermal optimum, and a unique ID.  
"""
abstract type AbstractSpecies end

struct Species <: AbstractSpecies
    np::Float64         # niche position
    r::Float64          # feeding range
    c_fr::Float64       # center of the feeding range
    id::UUID            # unique species ID
end

abstract type AbstractCommunity end

"""
    Community

Type containing parameters defining a community including its interaction matrix and a vector of species within.
"""
struct Community <: AbstractCommunity
    N::Int                      # Number of species
    A::Matrix{Float64}          # Interaction matrix
    sp::Vector{Species}         # species vector (includes the Species objects for each species)
    ids::Vector{UUID}           # unique species IDs
    np::Vector{Float64}         # niche parameter
    R::Float64
end

Base.show(io::IO, com::AbstractCommunity) =
    print(io, typeof(com), " N:", length(com.sp))

# --------------------------------------------------------------------------------------------------------
# Constructors for Species
# --------------------------------------------------------------------------------------------------------

"""
    species()

Constructor function for Sp type. Randomly generates niche parameters given a connectance value

Input:
C   : connectance value

Optional:
np  : niche position (default: randomly assigned)

Output:
Species object with all parameters as a vector of Species
"""
function species(C::Float64; np::Float64 = rand())
    # niche parameters
    # beta controls the feeding range distribution
    β = (1 / (2*C)) - 1
    
    # range of prey the species eats
    # picked randomly from a Beta distribution and linked to the niche parameter
    r = np*rand(Beta(1.0, β))

    # center of feeding range
    c_fr = rand(Uniform(r/2, np))

    # add a unique species ID
    uuid = UUIDs.uuid1()
    return Species(np, r, c_fr, uuid)
end

"""
    species(C::Float64, N::Int64)

Generate multiple species at the same time (similar to the one species version)
based on connectance.

Input:
C   : connectance value
N   : (Integer) number of species 

Optional:
np  : niche position (default: randomly assigned)
"""
function species(C::Float64, N::Int64; np::Vector{Float64} = rand(N))
    # control parameter for the feeding range distribution calculated from connectance C
    β = (1 / (2*C)) - 1

    # feeding range for each species picked from a Beta distribution
    r = np .* rand(Beta(1.0, β), N)

    # center of feeding range for each species
    c_fr = rand.(Uniform.(r/2, n))

    # add a unique species ID for each species
    uuid = [UUIDs.uuid1() for i = 1:N]

    return Species.(np, r, c_fr, uuid)
end

species(C::Float64, N::Int64; n::Vector{Float64} = rand(N)) = species(C, rand(N), N, n = n)

# --------------------------------------------------------------------------------------------------------
# Communities
# --------------------------------------------------------------------------------------------------------

"""
    community(sp_vec::Vector{Species}; R::Float64 = 42.0)

Generates an adjacency matrix for a given set of species using the niche model. 

Input:
sp_vec : vector of Species objects

Optional:
R      : (default = 42.0)

"""
function community(sp_vec::Vector{Species}; R::Float64 = 42.0)

    # Get number of species
    N = length(sp_vec)

    # prepare adjacency matrix as zero matrix
    A = zeros(N, N)

    # get vector for species niche parameter
    np = [sp_i.np for sp_i in sp_vec]

    # Loop over species. 
    # Set A[i, j] = 1 if species j is inside the range of i
    for (i, sp_i) = enumerate(sp_vec)
        for (j, sp_j) = enumerate(sp_vec)
            # if j is within range of i
            if sp_j.np < (sp_i.c_fr + (sp_i.r/2)) && sp_j.np > (sp_i.c_fr - (sp_i.r/2))
                A[i, j] = 1
            end
        end # end loop over j
    end # end loop over i

    # Add IDs
    ids = [x.id for x = sp_vec]

    # Create community as object
    com = Community(N, A, sp_vec, ids, np, R)

    # Check structure and remove mutual predation and cannibalism
    check_web!(com)

    return com
end

# --------------------------------------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------------------------------------

"""
    check_web!(com)

Clean up for the foodweb.
Remove any double and cannibalistic links from a community web. 
Double links represent mutual predation. Keep only the species with the bigger niche
as predator.
"""
function check_web!(com::Community)
    # For readability extract the adjacency matrix A 
    A = com.A
    # and vector of species objects
    sp_vec = com.sp

    # Find all double links (represents mutual predation: i eats j and j eats i)
    # If A[i, j] = 1 that means species i eats j
    # Multiplication of the upper and (transpose of the) lower triangular is only one, if there is a one
    # at both positions, aka i eats j and j eats i
    double_links = findall(UpperTriangular(A) .* LowerTriangular(A)' .!= 0)

    # Loop over and clean double links.
    # Compare niche parameter np; and keep the consumer with the larger np
    for link = double_links
        i, j = link.I
        if sp_vec[i].np > sp_vec[j].np
            A[j, i] = 0     # remove link where j eats i
        else
            A[i, j] = 0     # else remove i eats j
        end
    end

    # Remove cannibalism. 
    # Cannibals have a 1 on the diagonal.
    A[diagind(A)] .= 0

end

function check_iso(com::Community)

    # Find components
    g = SimpleDiGraph(com.A)

    components = connected_components(g)

    n_components = length(components)

    return n_components
    
end # end function

end # end module

