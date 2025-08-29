module FunctionalMotifs

# Import needed modules
include("MotifFinder.jl")
include("GeneralizedModel.jl")
include("Reactivity.jl")
include("Niche.jl")

using LinearAlgebra,
.MotifFinder,
.GeneralizedModel,
.Reactivity,
.Niche

export MotifResults, func_react_motifs, react_motifs_niche

"""
    MotifResults

A struct to hold all results from functional motif analysis.
z     : reactivity value of the most reactive instance of a motif 
inds  : indices of the most reactive instance of a motif 
props : proportion of system reactivity
"""
struct MotifResults
    z::Dict{Symbol, Vector{Float64}}
    inds::Dict{Symbol, Vector{Vector{Any}}}
    props::Dict{Symbol, Vector{Float64}}
end

"""
    MotifResultsNiche

A struct to hold all results from functional motif analysis.
z     : reactivity value of the most reactive instance of a motif 
inds  : indices of the most reactive instance of a motif 
props : proportion of system reactivity
Inter : interaction matrix produced by the niche model
"""
struct MotifResultsNiche
    z::Dict{Symbol, Vector{Float64}}
    inds::Dict{Symbol, Vector{Vector{Any}}}
    props::Dict{Symbol, Vector{Float64}}
    Inter::Dict{Symbol, Vector{Matrix{Int64}}}
end

"""
    func_react_motifs(n, interactions)

    Run analysis for functional reactivity motifs on a given interaction matrix 'interactions' for 'n' iterations.
    Returns a 'MotifResults' object.

    Different motifs are addressed by keys.
    This includes: 
    H : symmetric Jacobian of the system
    sp : species
    pp : predator-prey
    mut : mutualism (inside the apparent competition motif)
    ac : apparent competition
    ec : exploitative competition
    omni : omnivory
    ch3 : food chain with three species
    ch4 : food chain with four species
    di : diamond (combines ac and ec motifs)
    c3ch : combined tri-trophic food chains that share the same top-predator
"""
function func_react_motifs(n::Int, interactions::Matrix{Int64})

# Safety first
    if !(eltype(interactions) == Int64 && ndims(interactions) > 1)
        throw(ArgumentError("Interactions must be a Matrix{Int64} with more than one dimension."))
    end

    # Define motif keys
    # This includes: 
    # H : symmetric Jacobian of the system
    # sp : species
    # pp : predator-prey
    # mut : mutualism (inside the apparent competition motif)
    # ac : apparent competition
    # ec : exploitative competition
    # omni : omnivory
    # ch3 : food chain with three species
    # ch4 : food chain with four species
    # di : diamond (combines ac and ec motifs)
    # c3ch : combined tri-trophic food chains that share the same top-predator
    motif_keys = [:H, :sp, :pp, :mut, :ac, :ec, :ch3, :omni, :ch4, :di, :c3ch]

    # Find indices of all motifs ------------------------------------------------
    ind_pp, _ = MotifFinder.find_predprey(interactions)
    prey = getindex.(ind_pp, 1)
    pred = getindex.(ind_pp, 2)
    inds_pp = hcat(prey, pred)

    inds_ch3, _ = MotifFinder.find_threechain(interactions)
    inds_ac, _, inds_om, _ = MotifFinder.find_appcomp(interactions)
    inds_ec, _, inds_om2, _ = MotifFinder.find_explcomp(interactions)
    inds_omni, _ = MotifFinder.find_omni(inds_om, inds_om2)
    inds_mut, _ = MotifFinder.find_mutualism(interactions)
    inds_di, _ = MotifFinder.find_diamond(interactions)
    inds_ch4, _ = MotifFinder.find_fourchain(interactions)
    inds_c3ch, _ = MotifFinder.find_comb3chains(interactions)

    # Initialize result storage with preallocation ------------------------------
    z = Dict{Symbol, Vector{Float64}}()             # reactivity values
    inds = Dict{Symbol, Vector{Vector{Any}}}()      # indices of most reactive motif
    props = Dict{Symbol, Vector{Float64}}()         # proportion of system reactivity

    for k in motif_keys
        z[k] = Vector{Float64}(undef, n)
    end

    # exclude :H from inds and props as these only include values representing motifs, not the entire system
    for k in motif_keys[2:end]  
        inds[k] = Vector{Vector{Any}}(undef, n)
    end

    for k in motif_keys[2:end] 
        props[k] = Vector{Float64}(undef, n)
    end

    # get matrix size
    N = size(interactions)[1]

    for i in 1:n

        # Setup symmetric Jacobian
        H = zeros(N, N)

        # make sure only reactive systems are used
        # system is reactive if the largest eigenvalue of the symmetric matrix is positive
        while maximum(eigvals(H)) <= 0
            # Generate Jacobian with a Generalized Foodweb Model (the function only generates stable systems)
            J = GeneralizedModel.gen_J(interactions)
            # Calculate the symmetric Jacobian
            H = (J + transpose(J)) / 2
        end

        # System reactivity
        reac_H, _ = Reactivity.calc_reactivity(H)                       # calculates all eigenvalues, we only need the leading one
        z[:H][i] = reac_H

        # Find the most reactive species and save results
        max_sp, idx_sp = findmax(diag(H))
        z[:sp][i] = max_sp
        inds[:sp][i] = [idx_sp]
        props[:sp][i] = z[:sp][i] / z[:H][i]

        # Motif-specific reactivities
        reac_pp, idx_pp = Reactivity.calc_motif_reactivity(H, inds_pp) # calculates the reactivity of all motifs, we only need the most reactive
        z[:pp][i] = maximum(reac_pp)
        inds[:pp][i] = idx_pp
        props[:pp][i] = z[:pp][i] / z[:H][i]

        reac_mut, idx_mut = Reactivity.calc_motif_reactivity(H, inds_mut)
        z[:mut][i] = maximum(reac_mut)
        inds[:mut][i] = idx_mut
        props[:mut][i] = z[:mut][i] / z[:H][i]

        reac_ac, idx_ac = Reactivity.calc_motif_reactivity(H, inds_ac)
        z[:ac][i] = maximum(reac_ac)
        inds[:ac][i] = idx_ac
        props[:ac][i] = z[:ac][i] / z[:H][i]

        reac_ec, idx_ec = Reactivity.calc_motif_reactivity(H, inds_ec)
        z[:ec][i] = maximum(reac_ec)
        inds[:ec][i] = idx_ec
        props[:ec][i] = z[:ec][i] / z[:H][i]

        reac_ch3, idx_ch3 = Reactivity.calc_motif_reactivity(H, inds_ch3)
        z[:ch3][i] = maximum(reac_ch3)
        inds[:ch3][i] = idx_ch3
        props[:ch3][i] = z[:ch3][i] / z[:H][i]

        reac_omni, idx_omni = Reactivity.calc_motif_reactivity(H, inds_omni)
        z[:omni][i] = maximum(reac_omni)
        inds[:omni][i] = idx_omni
        props[:omni][i] = z[:omni][i] / z[:H][i]

        reac_ch4, idx_ch4 = Reactivity.calc_motif_reactivity(H, inds_ch4)
        z[:ch4][i] = maximum(reac_ch4)
        inds[:ch4][i] = idx_ch4
        props[:ch4][i] = z[:ch4][i] / z[:H][i]

        reac_di, idx_di = Reactivity.calc_motif_reactivity(H, inds_di)
        z[:di][i] = maximum(reac_di)
        inds[:di][i] = idx_di
        props[:di][i] = z[:di][i] / z[:H][i]

        reac_c3ch, idx_c3ch = Reactivity.calc_motif_reactivity(H, inds_c3ch)
        z[:c3ch][i] = maximum(reac_c3ch)
        inds[:c3ch][i] = idx_c3ch
        props[:c3ch][i] = z[:c3ch][i] / z[:H][i]

        # See if everything is running properly, give a little update every 1000 iterations
        if i % 1000 == 0
            println("Iteration: $i")
        end
    end

    return MotifResults(z, inds, props)
end # end function

"""
    react_motifs_niche(n, interactions)

    Run analysis for functional reactivity motifs on a given interaction matrix 'interactions' for 'n' iterations.
    Returns a 'MotifResults' object.

    Input:
    n : set size (number of iterations)
    N : number of species

    Optional:
    C : connectance (default 0.2)

    Different motifs are addressed by keys.
    This includes: 
    H : symmetric Jacobian of the system
    sp : species
    pp : predator-prey
    mut : mutualism (inside the apparent competition motif)
    ac : apparent competition
    ec : exploitative competition
    omni : omnivory
    ch3 : food chain with three species
    ch4 : food chain with four species
    di : diamond (combines ac and ec motifs)
    c3ch : combined tri-trophic food chains that share the same top-predator
"""
function react_motifs_niche(n::Int, N::Int; C::Float64 = 0.2)

    # Define motif keys
    # This includes: 
    # H : symmetric Jacobian of the system
    # sp : species
    # pp : predator-prey
    # mut : mutualism (inside the apparent competition motif)
    # ac : apparent competition
    # ec : exploitative competition
    # omni : omnivory
    # ch3 : food chain with three species
    # ch4 : food chain with four species
    # di : diamond (combines ac and ec motifs)
    # c3ch : combined tri-trophic food chains that share the same top-predator
    motif_keys = [:H, :sp, :pp, :mut, :ac, :ec, :ch3, :omni, :ch4, :di, :c3ch]

    # Initialize result storage with preallocation ------------------------------
    z = Dict{Symbol, Vector{Float64}}()             # reactivity values
    inds = Dict{Symbol, Vector{Vector{Any}}}()      # indices of most reactive motif
    props = Dict{Symbol, Vector{Float64}}()         # proportion of system reactivity
    Inter = Dict{Symbol, Vector{Matrix{Int64}}}() 

    Inter[:H] = Vector{Matrix{Int64}}(undef, n)

    for k in motif_keys
        z[k] = Vector{Float64}(undef, n)
    end

        # exclude :H from inds and props as these only include values representing motifs, not the entire system
    for k in motif_keys[2:end]  
        inds[k] = Vector{Vector{Any}}(undef, n)
    end

    for k in motif_keys[2:end] 
        props[k] = Vector{Float64}(undef, n)
    end

    i = 1

    # Create n communities
    while i < n + 1

        # Create species vector based on a niche model approach
        # In short: species are ordered based on a niche parameter,
        # a center and range for the feeding range are picked for each species
        # species gets a link to all species inside the range
        # mutual predation is removed - the species with the larger niche value is the consumer
        # cannibalism is removed
        sp_vec = [Niche.species(C) for i = 1:N]

        # Combine to community/network
        com = Niche.community(sp_vec)

        n_components = 42

        while n_components > 1
            sp_vec = [Niche.species(C) for i = 1:N]

            # Combine to community/network
            com = Niche.community(sp_vec)

            n_components = Niche.check_iso(com)
        end

        # Convert the interaction matrix
        A = convert(Matrix{Int64}, com.A)
        Inter[:H][i] = A

        # Create empty Hermetian
        H = zeros(size(A))
        # Get leading eigenvalue
        ev = maximum(eigvals(H))

        # Only create reactive systems
        while ev <= 0   
            # Generate a Jocobian for this network structure
            J = GeneralizedModel.gen_J(A)

            # Calculate the Hermetian matrix
            H = (J + transpose(J))/2;

            # Get leading eigenvalue
            ev = maximum(eigvals(H))
        end

        # Find indices of all motifs ------------------------------------------------
        ind_pp, _ = MotifFinder.find_predprey(A)
        prey = getindex.(ind_pp, 1)
        pred = getindex.(ind_pp, 2)
        inds_pp = hcat(prey, pred)

        inds_ch3, n_ch3 = MotifFinder.find_threechain(A)
        inds_ac, n_ac, inds_om, _ = MotifFinder.find_appcomp(A)
        inds_ec, n_ec, inds_om2, _ = MotifFinder.find_explcomp(A)
        inds_omni, n_omni = MotifFinder.find_omni(inds_om, inds_om2)
        inds_mut, _ = MotifFinder.find_mutualism(A)
        inds_di, n_di = MotifFinder.find_diamond(A)
        inds_ch4, n_ch4 = MotifFinder.find_fourchain(A)
        inds_c3ch, n_c3ch = MotifFinder.find_comb3chains(A)
        
        # Reactivities --------------------------------------------------------------
        # System reactivity
        reac_H, _ = Reactivity.calc_reactivity(H)                       # calculates all eigenvalues, we only need the leading one
        z[:H][i] = reac_H

        # Find the most reactive species and save results
        max_sp, idx_sp = findmax(diag(H))
        z[:sp][i] = max_sp
        inds[:sp][i] = [idx_sp]
        props[:sp][i] = z[:sp][i] / z[:H][i]

        # Motif-specific reactivities
        reac_pp, idx_pp = Reactivity.calc_motif_reactivity(H, inds_pp) # calculates the reactivity of all motifs, we only need the most reactive
        z[:pp][i] = maximum(reac_pp)
        inds[:pp][i] = idx_pp
        props[:pp][i] = z[:pp][i] / z[:H][i]

        # To account for the possibility that a motif is not part of the network, only safe them if they are.
        # If a motif is not part of a certain community, replace values by -999 for easy filtering
        if n_ac != 0
            reac_mut, idx_mut = Reactivity.calc_motif_reactivity(H, inds_mut)
            z[:mut][i] = maximum(reac_mut)
            inds[:mut][i] = [idx_mut]
            props[:mut][i] = z[:mut][i] / z[:H][i]
        else
            z[:mut][i] = -999
            inds[:mut][i] = [-999, -999]
            props[:mut][i] = -999
        end

        if n_ac != 0
            reac_ac, idx_ac = Reactivity.calc_motif_reactivity(H, inds_ac)
            z[:ac][i] = maximum(reac_ac)
            inds[:ac][i] = [idx_ac]
            props[:ac][i] = z[:ac][i] / z[:H][i]
        else
            z[:ac][i] = -999
            inds[:ac][i] = [-999, -999, -999]
            props[:ac][i] = -999
        end

        if n_ec != 0
            reac_ec, idx_ec = Reactivity.calc_motif_reactivity(H, inds_ec)
            z[:ec][i] = maximum(reac_ec)
            inds[:ec][i] = [idx_ec]
            props[:ec][i] = z[:ec][i] / z[:H][i]
        else
            z[:ec][i] = -999
            inds[:ec][i] = [-999, -999, -999]
            props[:ec][i] = -999
        end 

        if n_ch3 != 0
            reac_ch3, idx_ch3 = Reactivity.calc_motif_reactivity(H, inds_ch3)
            z[:ch3][i] = maximum(reac_ch3)
            inds[:ch3][i] = [idx_ch3]
            props[:ch3][i] = z[:ch3][i] / z[:H][i]
        else
            z[:ch3][i] = -999
            inds[:ch3][i] = [-999, -999, -999]
            props[:ch3][i] = -999
        end
    
        if n_omni != 0
            reac_omni, idx_omni = Reactivity.calc_motif_reactivity(H, inds_omni)
            z[:omni][i] = maximum(reac_omni)
            inds[:omni][i] = [idx_omni]
            props[:omni][i] = z[:omni][i] / z[:H][i]
        else
            z[:omni][i] = -999
            inds[:omni][i] = [-999, -999, -999]
            props[:omni][i] = -999
        end

        if n_ch4 != 0
            reac_ch4, idx_ch4 = Reactivity.calc_motif_reactivity(H, inds_ch4)
            z[:ch4][i] = maximum(reac_ch4)
            inds[:ch4][i] = [idx_ch4]
            props[:ch4][i] = z[:ch4][i] / z[:H][i]
        else
            z[:ch4][i] = -999
            inds[:ch4][i] = [-999, -999, -999, -999]
            props[:ch4][i] = -999
        end

        if n_di != 0
            reac_di, idx_di = Reactivity.calc_motif_reactivity(H, inds_di)
            z[:di][i] = maximum(reac_di)
            inds[:di][i] = [idx_di]
            props[:di][i] = z[:di][i] / z[:H][i]
        else
            z[:di][i] = -999
            inds[:di][i] = [-999, -999, -999, -999]
            props[:di][i] = -999
        end

        if n_c3ch != 0
            reac_c3ch, idx_c3ch = Reactivity.calc_motif_reactivity(H, inds_c3ch)
            z[:c3ch][i] = maximum(reac_c3ch)
            inds[:c3ch][i] = [idx_c3ch]
            props[:c3ch][i] = z[:c3ch][i] / z[:H][i]
        else
            z[:c3ch][i] = -999
            inds[:c3ch][i] = [-999, -999, -999, -999, -999]
            props[:c3ch][i] = -999
        end

        i += 1

        # See if everything is running properly, give a little update every 1000 iterations
        if i % 1000 == 0
            println("Iteration: $i")
        end

    end

    return MotifResultsNiche(z, inds, props, Inter)
    
end

end # module