module MotifFinder

using Combinatorics
using DataFrames
using Distributions
using IterTools
using LinearAlgebra
using Random
using StatsBase

export find_predprey, find_mutualism, find_threechain, find_appcomp, find_explcomp, find_diamond, find_fourchain, find_bifan, find_comb3chains

# -------------------------------------------------------------------------------------------------------------------------------
# Functions included in this module are able to find and count given motifs in an interaction matrix of a food web.

# The interaction matrix must contain of ones and zeros only. Columns represent predation, rows stand for prey species. 
# A one indicates an interaction where the predator (column) species feeds on the prey (row) species.

# Available motifs: Predator-Prey, Mutualism, Three Chain, Apparent Competition, Exploitative Competition, Four Chain, Diamond
# Included functions: 
#                   find_predprey() 
#                   find_mutualism()
#                   find_threechain()
#                   find_appcomp()
#                   find_explcomp()
#                   find_fourchain()
#                   find_diamond()
#                   find_bifan()
#                   find_comb3chains()

# General Function Structure:
# Function call: 
# ind, n = find_...(input)

# Input:            
#                   interactions: matrix classifying the interactions (works with all functions)
#                   or 
#                   inds:         vector of CartesianIndices (in motifs with > 2 species)
#                   or:
#                   inds:         matrix with the indices of a three species motif (specifically for mutualism)
# Output:
#                   ind: matrix with indices of species in the motif
#                   n:   number of found motifs in the food web

# NOTE: The function find_mutualism() is primarily used to find indirect mutualistic pairs in the apparent competition motif!
# -------------------------------------------------------------------------------------------------------------------------------

"""
    find_predprey()

    Find and count predator-prey pairs. Returns indices and the number of found motifs

    Input: 
    interactions :
    Interaction matrix (Matrix{Int}) with 0s and 1s that mark who eats whom

    Output: 
    ind_predprey : CartesianIndices of all predator-prey pairs (first index is the prey species, second the predator)
    n_predprey   : number of predator-prey pairs
"""
function find_predprey(interactions::Matrix{Int})

    # Find all predator-prey pairs
    ind_predprey = findall(interactions .== 1);

    # Count motifs
    n_predprey = length(ind_predprey);

    return ind_predprey, n_predprey

end # end function

"""
    find_mutualism()

    Find and count mutualistic pairs in apparent competition motifs

    Input: 
    Interaction matrix or CartesianIndices with all predator-prey pairs

    Output: 
    ind_mutual : matrix with indices of mutualistic pairs
    n_mutual   : number of mutualisitc pairs
"""
function find_mutualism(input)
    if isa(input, Matrix) && size(input, 1) > 3
        # Input is an interaction matrix. First calculate all predator-prey pairs.
        ind_appcomp, _ = find_appcomp(input);
    else
        ind_appcomp = input;
    end

    # Get the first two columns with the pair of preys in the motif
    pairs = ind_appcomp[:, 1:2];

    # Only use unique rows to assure that no pair is counted multiple times
    ind_mutual = unique(pairs, dims = 1);

    # Count motifs
    n_mutual = size(ind_mutual, 1);
    
    return ind_mutual, n_mutual
end # end function

"""
    find_threechain()

    Find the three chains/tri-trophic food chain (food chains of three species) in a food web.

    Input: 
    Interaction matrix or CartesianIndices with all predator-prey pairs

    Output:
    ind_threechain : matrix with indices of species, 1 column: lowest trophic level, 2 column : intermediate level, 3 column: top level
    n_threechain   : number of three chains (or tri-trophic food chains) 
"""
function find_threechain(input)
    if isa(input, Matrix)
        # Input is an interaction matrix. First calculate all predator-prey pairs.
        ind_predprey, _ = find_predprey(input);
    else
        ind_predprey = input;
    end

    # If indices of all predator-prey pairs are available, create sets of indices for species in a three chain.
    # First entry in the CartesianIndices is the prey index, second is the predator index. 
    prey = getindex.(ind_predprey, 1);
    pred = getindex.(ind_predprey, 2);

    pp = hcat(prey, pred);

    # Compare both to check if a species is both a predator and a prey to some species.
    common_entries = collect(intersect(Set(pred), Set(prey)));
    
    # Set up empty matrix for three chain indices
    ind_threechain = zeros(Int, 0, 3);

    for n in common_entries
        # Find indices where n is the prey species
        ind_prey = findall(x -> x == n, prey);

        # Find indices where n is the predator species
        ind_pred = findall(x -> x == n, pred);

        # Get predator-prey pair where n is the predator
        pp_prey = ind_predprey[ind_prey];
        pp_pred = ind_predprey[ind_pred];

        # Extract elements from each index and construct the result matrix using array comprehension
        inds = hcat([[pp_pred_elem[1], pp_pred_elem[2], pp_prey_elem[2]] for pp_pred_elem in pp_pred, pp_prey_elem in pp_prey]...);

        inds_unique = unique(inds, dims = 2);

        ind_threechain = [ind_threechain; inds_unique'];
    end

    # Remove rows where prey and top predator are a predator-prey pair
    # Identify rows in mat1 to be removed
    rows_to_remove = [row_in_matrix(ind_threechain[i, [1, 3]], pp) for i in 1:size(ind_threechain, 1)]

    # Filter the rows to be kept
    ind_threechain = ind_threechain[.!rows_to_remove, :]

    # Remove potential duplicates
    ind_threechain = rem_rows_with_dupl(ind_threechain);

    # Count motifs
    n_threechain = size(ind_threechain, 1);

    return ind_threechain, n_threechain

end # end function

"""
    find_explcomp()

    Find the exploitative competition motif where two predator species share the same prey.

    Input: 
    Interaction matrix or CartesianIndices with all predator-prey pairs

    Output:
    ind_explcomp : matrix with indices of species in the motif, first two columns are the prey species
    ind_omni     : matrix with indices of omnivory motifs
    n_explcomp   : number of exploitative competition motifs
    n_omni       : number of omnivory motifs
"""
function find_explcomp(input)
    if isa(input, Matrix)
        # Input is an interaction matrix. First calculate all predator-prey pairs.
        ind_predprey, _ = find_predprey(input);
    else
        ind_predprey = input;
    end

    # If indices of all predator-prey pairs are available, create sets of indices for species in a three chain.
    # First entry in the CartesianIndices is the prey index, second is the predator index. 
    prey = getindex.(ind_predprey, 1);
    pred = getindex.(ind_predprey, 2);

    # Combine to matrix
    predprey = hcat(prey, pred);

    # Create a set
    predprey_set = Set(sort.(eachrow(predprey)));

    # Count how often the species appear. Species that appear more than once are part of the motif
    count_preds = countmap(pred);

    # Isolate the counts
    counts = [count_preds[key] for key in collect(keys(count_preds))];

    # Get the species that appears more than once
    species = collect(keys(count_preds))[findall(x -> x > 1, counts)];

    # Set up empty matrix for three chain indices
    ind_explcomp = zeros(Int, 0, 3);
    ind_omni = zeros(Int, 0, 3);

    for n in species
        # Get the indices of all predators n
        ind_pred = findall(x -> x == n, pred);

        # Get the preys with the same predator
        preys = prey[ind_pred];

        # Generate all possible pairs of combinations with two different prey species and add the predator as third column
        combis = hcat([[preys[i], preys[j], n] for i in eachindex(preys), j in eachindex(preys) if i < j]...);

        # Transpose
        combis = combis';
        
        uniq_combis = unique(combis, dims = 1);

        # Get the two predator species
        dummy = uniq_combis[:, 1:2];

        # Check for rows that don't already appear in the predator-prey pairs
        keep = [row ∉ predprey_set for row in eachrow(sort(dummy, dims = 2))];

        # Combinations where the preys don't form a predator-prey pair are counted as apparent competition
        ind_explcomp = [ind_explcomp; combis[keep, :]];
        # The others are the omnivory motif
        ind_omni = [ind_omni; combis[.!keep, :]];
        
    end

    ind_omni = rem_rows_with_dupl(ind_omni);
    ind_explcomp = rem_rows_with_dupl(ind_explcomp);
 
    # Count motifs
    n_explcomp = size(ind_explcomp, 1);
    n_omni = size(ind_omni, 1);

    return ind_explcomp, n_explcomp, ind_omni, n_omni
end # end function

"""
    find_appcomp()

    Find the apparent competition motif where two prey species share the same predator.

    Input: 
    Interaction matrix or CartesianIndices with all predator-prey pairs

    Output:
    ind_appcomp : matrix with the indices of the apparent competition motif, first two columns are the predator species
    ind_omni    : matrix with indices of the omnivory motif
    n_appcomp   : number of apparent competition motifs
    n_omni      : number of apparent competition motifs
"""
function find_appcomp(input)
    if isa(input, Matrix)
        # Input is an interaction matrix. First calculate all predator-prey pairs.
        ind_predprey, _ = find_predprey(input);
    else
        ind_predprey = input;
    end

    # If indices of all predator-prey pairs are available, create sets of indices for species in a three chain.
    # First entry in the CartesianIndices is the prey index, second is the predator index. 
    prey = getindex.(ind_predprey, 1);
    pred = getindex.(ind_predprey, 2);

    # Combine to matrix
    predprey = hcat(prey, pred);

    # Create a set
    predprey_set = Set(sort.(eachrow(predprey)));

    # Count how often the species appear. Species that appear more than once are part of the motif
    count_preys = countmap(prey);

    # Isolate the counts
    counts = [count_preys[key] for key in collect(keys(count_preys))];

    # Get the species that appears more than once
    species = collect(keys(count_preys))[findall(x -> x > 1, counts)];

    # Set up empty matrix for three chain indices
    ind_appcomp = zeros(Int, 0, 3);
    ind_omni = zeros(Int, 0, 3);

    for n in species
        # Get the indices of all preys n
        ind_prey = findall(x -> x == n, prey);

        # Get the preys with the same predator
        preds = pred[ind_prey];

        # Generate all possible pairs of combinations with two different prey species and add the predator as third column
        combis = hcat([[preds[i], preds[j], n] for i in eachindex(preds), j in eachindex(preds) if i < j]...);

         # Transpose
        combis = combis';

        # Get the two prey species
        dummy = combis[:, 1:2];

        # Check for rows in mat1 not present in mat2
        not_in_matrix = [row ∉ predprey_set for row in eachrow(sort(dummy, dims = 2))];

        ind_appcomp = [ind_appcomp; combis[not_in_matrix, :]];

        ind_omni = [ind_omni; combis[.!not_in_matrix, :]];

    end
    # Remove potential duplicates
    ind_omni = rem_rows_with_dupl(ind_omni);
    ind_appcomp = rem_rows_with_dupl(ind_appcomp);
 
    # Count motifs
    n_appcomp = size(ind_appcomp, 1);
    n_omni = size(ind_omni, 1); 

    return ind_appcomp, n_appcomp, ind_omni, n_omni
end # end function


"""
    find_omni()

    Find the omnivory motifs of three species in a food web.

    Input: Interaction matrix or CartesianIndices with all predator-prey pairs

    Output:
    ind_omni : matrix with indices of species, 1 column: intermediate level, 2 colum : top level, 3 column: lowest level
    n_omni   : number of three chains (or tri-trophic food chains) 
"""
function find_omni(args...)
    omni_ac = args[1];
    omni_ec = args[2];

    omni = vcat(omni_ac, omni_ec);

    ind_omni = remove_duplicate_and_permuted_rows(omni)

    n_omni = size(ind_omni, 1); 

    return ind_omni, n_omni
    
end

"""
    find_diamond()

    Find the diamond motif. Diamond motifs are a combination of apparent competition and exploitative competition. 
    The prey in the apparent competition are the same as the predators in the exploitative competition motif.

    Output:
    ind_diamond : matrix with the indices of the diamond motif, first two columns are the middle species, third is the prey
                  of the exploitative competition motif and the fourth the top predator in this motif
    n_diamond   : number of diamond motifs
"""
function find_diamond(args...)
    if length(args) == 1
        # Input is an interaction matrix. First calculate all predator-prey pairs.
        ind_appcomp, _ = find_appcomp(args[1]);
        ind_explcomp, _ = find_explcomp(args[1]);
    elseif length(args) == 2
        ind_appcomp = args[1];
        ind_explcomp = args[2];
    else
        println("Function does not support this form of input. Please read the module introduction to this function.")
    end

    # Take the first two columns of both matrices and compare the rows
    preys = ind_appcomp;
    preds = ind_explcomp;

    # Define the columns that should be compared
    cols = (1, 2);

    # The diamond is a combination of apparent competition and exploitative competition.
    # Merge the two motifs depending on the first two columns.
    ind_diamond = [[row1[1], row1[2], row1[3], row2[3]] for row1 in eachrow(preys), 
                    row2 in eachrow(preds) if row1[cols[1]] == row2[cols[1]] && row1[cols[2]] == row2[cols[2]]];

    ind_diamond = hcat(ind_diamond...)';

    ind_diamond = Matrix(ind_diamond);

    # Remove potential duplicates
    ind_diamond = rem_rows_with_dupl(ind_diamond);

    # Count motifs
    if size(ind_diamond, 2) == 0
        n_diamond = 0;
    else
        n_diamond = size(ind_diamond, 1);
    end
    
    return ind_diamond, n_diamond
end # end function

"""
    find_fourchain()

    Find and count food chains with four species

    Input: 
    Interaction matrix or CartesianIndices with all predator-prey pairs

    Output: 
    ind_fourchain : matrix with indices of four-species-food-chain motifs (sorted from first column = lowest level prey to 
                    last column = top predator in the motif)
    n_fourchain   : number of four-species-food-chain motifs
"""
# Function to find and count all four chains
function find_fourchain(args...)
    if length(args) == 1
        # Input is an interaction matrix. First calculate all predator-prey pairs.
        ind_pp, _ = find_predprey(args[1]);
        ind_three, _ = find_threechain(args[1]);
    elseif length(args) == 2
        ind_pp = args[1];
        ind_three = args[2];
    else
        println("Function does not support this form of input. Please read the module introduction to this function.")
    end

    # Use the three chains and expand by a fourth species
    # Either the species on the lowest trophic level (first column in the matrix) is predator to another species 
    # or the species on the highest level (third column) is prey to another species
    lowestlvl = ind_three[:, 1];
    highestlvl = ind_three[:, 3];

    # Convert pred-prey pairs to matrix
    prey = getindex.(ind_pp, 1);
    pred = getindex.(ind_pp, 2);

    matrix = hcat(prey, pred);

    # Initialise result matrix
    ind_fourchain = zeros(Int, 0, 4);

    for n in lowestlvl
        ind = findall(x -> x == n, pred);

        dummy = matrix[ind, :];

        # Define the columns that should be compared
        cols = (1, 2);

        # Merge the predator-prey pairs to a four chain
        ind_four = [[row1[1], row2[1], row2[2], row2[3]] for row1 in eachrow(dummy), 
                        row2 in eachrow(ind_three) if row1[cols[2]] == row2[cols[1]]];

        if !isempty(ind_four)
            ind_fourchain = [ind_fourchain; hcat(ind_four...)'];
        end

    end

    # Repeat for the highest level
    for n in highestlvl
        ind = findall(x -> x == n, prey);

        dummy = matrix[ind, :];

        # Define the columns that should be compared
        cols = (1, 3);

        ind_four = [[row2[1], row2[2], row2[3], row1[2]] for row1 in eachrow(dummy), 
                        row2 in eachrow(ind_three) if row1[cols[1]] == row2[cols[2]]];

        if !isempty(ind_four)
            ind_fourchain = [ind_fourchain; hcat(ind_four...)'];
        end

    end

    # Get rid off double counts
    ind_fourchain = unique(ind_fourchain, dims = 1);

    result = zeros(Int, 0, 4);

    for i in eachrow(ind_fourchain)

        v = vec(i)
        if length(unique(v)) == 4
            result = [result; i'];
        end

    end

    ind_fourchain = result;

            # Identify rows in mat1 to be removed
            rows_to_remove = [row_in_matrix(ind_fourchain[i, [1,4]], matrix) for i in 1:size(ind_fourchain, 1)]
            ind_fourchain = ind_fourchain[.!rows_to_remove, :]

                        # Identify rows in mat1 to be removed
                        rows_to_remove = [row_in_matrix(ind_fourchain[i, [2,4]], matrix) for i in 1:size(ind_fourchain, 1)]
                        ind_fourchain = ind_fourchain[.!rows_to_remove, :]
            
    # Count the motifs
    n_fourchain = size(ind_fourchain, 1);

    return ind_fourchain, n_fourchain
    
end # end function

"""
    find_bifan()

    Find and count bifan motifs with four species. A bifan is a combination of two exploitative competition motifs where
    two prey species share the same two predators.

    Input: 
    Interaction matrix or CartesianIndices with all predator-prey pairs

    Output: 
    ind_bifan : matrix with indices of bifan motifs, first two columns are the prey species, third and fourth column are the predators
    n_bifan   : number of bifan motifs
"""
function find_bifan(args...)
    if length(args) == 2 && size(args[1], 2) > 3
        # Input is an interaction matrix. First get all exploitative competition sets.
        ind_ec, _ = find_explcomp(args[1]);
        ind_pp = args[2];

    elseif length(args) == 2
        ind_ec = args[1];
        ind_pp = args[2];
    else
        println("Function does not support this form of input. Please use either an interaction matrix or a nx3 matrix of indices!")
    end

    prey = getindex.(ind_pp, 1);
    pred = getindex.(ind_pp, 2);

    ind_pp = hcat(prey, pred);

    # The bi-fan motif is a combination of two exploitative competition motifs that share the same predators
    # Define the columns to compare
    cols = (1, 2);

    # Generate bifan motifs
    bifan = [[row1[1], row1[2], row1[3], row2[3]] for row1 in eachrow(ind_ec),
            row2 in eachrow(ind_ec) if row1[cols[1]] == row2[cols[1]] && row1[cols[2]] == row2[cols[2]] && row1[3] != row2[3]];

    # Convert to a matrix if not empty
    bifan = length(bifan) > 0 ? hcat(bifan...)' : [];

    # Convert rows of the predator-prey indices into a set of sorted tuples for fast lookup
    rows_pp = Set(Tuple.(sort.(eachrow(ind_pp))));

    # Identify rows in bifan to keep (considering only the last two columns)
    keep = [i for (i, row) in enumerate(eachrow(bifan)) if Tuple(sort(row[end-1:end])) ∉ rows_pp];

    # Just keep the rows where there is no predator-prey pair
    bifan = length(keep) > 0 ? bifan[keep, :] : [];

    # Remove those where the two preys are a predator-prey pair
    # Identify rows in mat1 to be removed
    rows_to_remove = [row_in_matrix(bifan[i, 1:2], ind_pp) for i in 1:size(bifan, 1)]

    # Filter the rows to be kept
    bifan = bifan[.!rows_to_remove, :]

    rows_to_remove = [row_in_matrix(bifan[i, [2, 1]], ind_pp) for i in 1:size(bifan, 1)]

    # Filter the rows to be kept
    bifan = bifan[.!rows_to_remove, :]

        # Identify rows in mat1 to be removed
        rows_to_remove = [row_in_matrix(bifan[i, 3:4], ind_pp) for i in 1:size(bifan, 1)]
        bifan = bifan[.!rows_to_remove, :]

    # Concatenate the sorted first two and last two columns of each vector
    sorted_bifan = [vcat(sort(v[1:2]), sort(v[3:4])) for v in eachrow(bifan)];

    # Only keep unique combinations
    ind_bifan = unique(sorted_bifan);

    # Convert to matrix if not empty
    ind_bifan = length(ind_bifan) > 0 ? hcat(ind_bifan...)' : [];

    # Count the motifs
    n_bifan = size(ind_bifan, 1);

    return ind_bifan, n_bifan
    
end # end function

"""
    find_comb3chains()

    Find and count combined three chain motifs with five species. Combined three chains are motifs where two three chain motifs
    share the same predator on the top level

    Input: 
    Interaction matrix or CartesianIndices with all predator-prey pairs

    Output: 
    ind_combtc : matrix with indices of combined three chain motifs, columns 1 to three give the first three chain in the order
                 top predator -> intermediate level -> lowest level, followed by the intermediate level and lowest level of the 
                 second three chain
    n_combtc   : number of bifan motifs
"""
function find_comb3chains(args...)
    if length(args) == 1 
        # Input is an interaction matrix. First calculate all predator-prey pairs and three chains.
        ind_pp, _ = MotifFinder.find_predprey(args[1]);
        ind_tc, _ = MotifFinder.find_threechain(ind_pp);
    elseif length(args) == 2
        ind_pp = args[1];
        ind_tc = args[2];
    else
        println("Function does not support this form of input. Please use either an interaction matrix or two matrices with indices for predator prey pairs and three chains.")
    end
    
    # Create matrix of predator prey pairs
    pp_pairs = hcat(getindex.(ind_pp, 1), getindex.(ind_pp, 2));

    # First two columns are the lower two trophic levels, third is the top predator in this motif
    low = ind_tc[:, 1];
    mid = ind_tc[:, 2];
    top = ind_tc[:, 3];

    # Get all top predators
    tops = unique(top)

    # Set up empty matrix for three chain indices
    ind_c3c = zeros(Int, 0, 5);
    
    for n in tops
        # Find indices where n is the top predator
        ind_tops = findall(x -> x == n, top);

        # Create all combinations
        for i in ind_tops
            for j in ind_tops
                if i != j
                    ind_c3c = [ind_c3c; top[i] mid[i] low[i] mid[j] low[j]]
                end
            end
        end

    end

    # Filter the matrix
    inds_c3c = remove_duplicate_and_permuted_rows(ind_c3c)

    # Remove those where the middle species are a predator-prey pair
    # Remove rows where prey and top predator are a predator-prey pair
    # Identify rows in mat1 to be removed
    rows_to_remove = [row_in_matrix(inds_c3c[i, [2, 4]], pp_pairs) for i in 1:size(inds_c3c, 1)]

    # Filter the rows to be kept
    inds_c3c = inds_c3c[.!rows_to_remove, :]

    rows_to_remove = [row_in_matrix(inds_c3c[i, [2, 5]], pp_pairs) for i in 1:size(inds_c3c, 1)]

    # Filter the rows to be kept
    inds_c3c = inds_c3c[.!rows_to_remove, :]

    rows_to_remove = [row_in_matrix(inds_c3c[i, [3, 4]], pp_pairs) for i in 1:size(inds_c3c, 1)]

    # Filter the rows to be kept
    inds_c3c = inds_c3c[.!rows_to_remove, :]

    rows_to_remove = [row_in_matrix(inds_c3c[i, [3, 5]], pp_pairs) for i in 1:size(inds_c3c, 1)]
    
    # Filter the rows to be kept
    inds_c3c = inds_c3c[.!rows_to_remove, :]

    # Remove double counts
    ind_combtc = rem_rows_with_dupl(inds_c3c);
     
    # Count motifs
    n_combtc = size(ind_combtc, 1);
    
    return ind_combtc, n_combtc

end

"""
    are_permutations()

    Function to check if two rows in a matrix are permutations of each other.

    Input: 
    row1 : first row to check
    row2 : second row to check

    Output:
    Boolean 1 if rows are permutations, 0 otherwise
"""
function are_permutations(row1, row2)
    return sort(row1) == sort(row2)
end

"""
    row_in_matrix()

    Checks if any row from a matrix appears as a permutation of a row of another matrix
"""
# Function to check if a row from mat1 appears as a permutation in mat2
function row_in_matrix(row, matrix)
    for i in eachrow(matrix)
        if are_permutations(row, i)
            return true
        end
    end
    return false
end

"""
    remove_duplicate_and_permuted_rows()

    Removes duplicates and permuted rows in a matrix (to avoid counting rows with the same species double)
"""
function remove_duplicate_and_permuted_rows(mat)
    nrows, ncols = size(mat)
    unique_rows = Matrix{eltype(mat)}(undef, 0, ncols)
    
    for i in 1:nrows
        row = mat[i, :]
        if !any(are_permutations(row, unique_rows[j, :]) for j in 1:size(unique_rows, 1))
            unique_rows = vcat(unique_rows, row')
        end
    end
    
    return unique_rows
end

function rem_rows_with_dupl(matrix::Matrix{Any})
    # Initialize an empty vector to store indices of rows to remove
    rows_to_remove = Int[]
    
    # Iterate through each row of the matrix
    for i in 1:size(matrix, 1)
        # Check if the row has any duplicate elements
        if length(unique(matrix[i, :])) < size(matrix, 2)
            push!(rows_to_remove, i)
        end
    end
    
    # Remove rows identified to be removed
    return matrix[setdiff(1:size(matrix, 1), rows_to_remove), :]
end 

function rem_rows_with_dupl(matrix::Matrix{Int64})
    # Initialize an empty vector to store indices of rows to remove
    rows_to_remove = Int[]
    
    # Iterate through each row of the matrix
    for i in 1:size(matrix, 1)
        # Check if the row has any duplicate elements
        if length(unique(matrix[i, :])) < size(matrix, 2)
            push!(rows_to_remove, i)
        end
    end
    
    # Remove rows identified to be removed
    return matrix[setdiff(1:size(matrix, 1), rows_to_remove), :]
end # end function

end # end module