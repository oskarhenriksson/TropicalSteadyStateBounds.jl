export minimal_presentation, 
    augmented_vertical_system

function augmented_vertical_system(C::QQMatrix, M::ZZMatrix, L::QQMatrix)

    KB, k, b = rational_function_field(QQ, "k"=>1:ncols(M), "b"=>1:nrows(L))
    R, x = polynomial_ring(KB, "x"=>1:nrows(M))

    steady_state = C*[k[j]*prod(x .^ M[:, j]) for j in 1:ncols(M)]
    conservation_laws = L*x .- b

    return [steady_state; conservation_laws]
end

"""
    minimal_presentation(C, M)

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

    julia> C_tilde, M_tilde = minimal_presentation(C, M);

    julia> C_tilde
    [   0           k[3]   -k[4]           k[5]]
    [k[1]   -k[2] - k[3]       0              0]
    [   0              0    k[4]   -k[5] - k[6]]

    julia> M_tilde
    [1   0   0   0]
    [0   0   1   0]
    [1   0   0   0]
    [0   0   1   0]
    [0   1   0   0]
    [0   0   0   1]

    ```
"""
function minimal_presentation(C::QQMatrix, M::ZZMatrix)
    columns = [M[:, i] for i in 1:ncols(M)]
    unique_columns = unique(columns)
    r = length(unique_columns)

    M_tilde = matrix(ZZ, hcat(unique_columns...))
    K, k = rational_function_field(QQ, "k"=>1:ncols(M))
    C_tilde = zero_matrix(K, nrows(C), r) 
    for i = 1:r
        indices = findall(c -> c == unique_columns[i], columns)
        for j in indices
            C_tilde[:, i] += k[j] .* C[:, j]
        end
    end
    return C_tilde, M_tilde
end