export monomial_reembedding

"""
    monomial_reembedding(C, M)

    Given a vertical system given by`C` and exponent matrix `M`, 
    compute the defining matrices for the monomial re-embedding.

    # Example
    ```jldoctest
    julia> C = matrix(QQ, [0 0 1 -1 1 0; 1 -1 -1 0 0 0; 0 0 0 1 -1 -1])
    [0    0    1   -1    1    0]
    [1   -1   -1    0    0    0]
    [0    0    0    1   -1   -1]

    julia> M = matrix(ZZ, [1 0 0 0 0 0; 0 0 0 1 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 1 1 0 0 0; 0 0 0 0 1 1])
    [1   0   0   0   0   0]
    [0   0   0   1   0   0]
    [1   0   0   0   0   0]
    [0   0   0   1   0   0]
    [0   1   1   0   0   0]
    [0   0   0   0   1   1]

    julia> C_tilde, M_tilde = monomial_reembedding(C, M);

    julia> C_tilde
    [   0           κ[3]   -κ[4]           κ[5]]
    [κ[1]   -κ[2] - κ[3]       0              0]
    [   0              0    κ[4]   -κ[5] - κ[6]]

    julia> M_tilde
    [1   0   0   0]
    [0   0   1   0]
    [1   0   0   0]
    [0   0   1   0]
    [0   1   0   0]
    [0   0   0   1]

    ```
"""
function monomial_reembedding(C::QQMatrix, M::ZZMatrix)
    columns = [M[:, i] for i in 1:ncols(M)]
    unique_columns = unique(columns)
    r = length(unique_columns)

    M_tilde = matrix(ZZ, hcat(unique_columns...))
    K, κ = rational_function_field(QQ, "κ"=>1:ncols(M))
    C_tilde = zero_matrix(K, nrows(C), r) 
    for i = 1:r
        indices = findall(c -> c == unique_columns[i], columns)
        for j in indices
            C_tilde[:, i] += κ[j] .* C[:, j]
        end
    end
    return C_tilde, M_tilde
end