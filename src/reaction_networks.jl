export augmented_vertical_system

"""
    augmented_vertical_system(rn)

    For a reaction network `rn`, compute the coefficient, exponent, and linear part matrices
    of the augmented vertically parametrized steady state system.

    # Example
    ```jldoctest
    julia> using Catalyst;
    
    julia> rn = @reaction_network begin
        k1, X1 --> X2
        k2, X2 --> X1
        k3, 2*X1 + X2 --> 3*X1
    end;

    julia> C, M, L = augmented_vertical_system(rn);

    julia> C
    [1   -1   -1]

    julia> M
    [1   0   2]
    [0   1   1]

    julia> L
    [1   1]
    ````

"""
function augmented_vertical_system(rn::ReactionSystem)

    @req all(ismassaction.(reactions(rn), Ref(rn))) "All reactions must have mass action kinetics"

    N = matrix(QQ, netstoichmat(rn)) #stoiciometric matrix
    M = matrix(ZZ, substoichmat(rn)) #kinetic matrix
    L = rref(matrix(QQ, conservationlaws(rn)))[2] #conserved quantities
    
    # Compute the matrix C (choose of linearly independent rows of N)
    first_nonzero_indices = [findfirst(!iszero, row) for row in eachrow(L)]
    C = N[setdiff(collect(1:nrows(N)), first_nonzero_indices), :] 

    return C, M, L
end



