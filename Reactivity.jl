module Reactivity

using LinearAlgebra

export calc_reactivity, calc_motif_reactivity

# --------------------------
# The module can be used to calculate the reactivity load z of a network and different kinds of motifs that are part of it.
# The reactivity load z is defined as the largest eigenvalue of the Hermitian
# --------------------------

"""
calc_reactivity()

Function to calculate the reactivity of a system based on the symmetric Jacobian matrix. The reactivity is defined as the 
largest eigenvalue of the matrix.

Input: 
S : symmetric Jacobian matrix of the system, Matrix{Float64}

Output:
z_S : reactivity
ind : index of the most reactive eigenvalue
"""
function calc_reactivity(S::Matrix{Float64})
    # Calculate the reactivity load z 
    z_S, ind = findmax(eigvals(S));

    return z_S, ind
end

"""
calc_motif_reactivity()

Function to calculate the reactivity of a system based on the symmetric Jacobian matrix. The reactivity is defined as the 
largest eigenvalue of the matrix.

Input: 
S    : symmetric Jacobian matrix of the system, Matrix{Float64}
inds : matrix with the indices of each occurrence of a certain motif

Output:
z_motif         : reactivity of each appearance of a certain motif 
inds_maxz_motif : indices of the nodes in the most reactive motif
"""
function calc_motif_reactivity(H::Matrix{Float64}, inds)

    # only calculate reactivity if there are any motifs of a certain kind
    try
        if inds === nothing
            return nothing
        end

        if isempty(inds) || is_empty_any_matrix(inds)
            return nothing
        else
        
        z_motif = zeros(1, size(inds, 1))
        k = 1

        if size(inds, 1) > 1
            for i in eachrow(inds)
                idx = sort(vec(i))
                z_motif[k] = maximum(eigvals(H[idx, idx]))
                k += 1
            end

            _, ind = findmax(z_motif)
            inds_maxz_motif = [inds[ind[2], :]]
        else
            idx = sort(vec(inds))
            H = reshape(H[idx, idx], size(inds, 2), size(inds, 2))
            z_motif[1] = maximum(eigvals(H))
            inds_maxz_motif = inds 
        end

        return z_motif, inds_maxz_motif

    end

    catch err
        return nothing
    end

end

""" 
is_empty_any_matrix()

Checks if a matrix of type Any is empty.
"""
function is_empty_any_matrix(mat)
    if !isa(mat, Matrix{Any})
        return false
    end
    
    nrows, ncols = size(mat)
    return ncols == 0 && nrows == 1
end

end # end module