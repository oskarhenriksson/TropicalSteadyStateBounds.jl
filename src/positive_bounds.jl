
export lower_bound_of_maximal_positive_steady_state_count,
    lower_bound_of_maximal_positive_root_count,
    lower_bound_of_maximal_positive_root_count_fixed_b_h


@doc raw"""
    lower_bound_of_maximal_positive_root_count_fixed_b_h(
    C::QQMatrix, M::ZZMatrix, L::QQMatrix,
    b_spec::Union{Vector{Int},Vector{QQFieldElem}}, 
    h::Union{Vector{Int},Vector{QQFieldElem}}; 
    TropB::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)

Compute a lower bound of the maximal positive root count for an augmented vertically parametrized system
given by the matrices `C`, `M` and `L`, given a mixed choice of constant terms `b_spec` and shift `h` of
the tropicalized binomial variety.

# Example
```jldoctest
julia> C = matrix(QQ, [1 -1 -1]);

julia> M = matrix(ZZ, [1 0 2; 0 1 1]);

julia> L = matrix(QQ, [1 1]);

julia> h = [37,97,18];

julia> b = [71];

julia> lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h)
3
```
"""
function lower_bound_of_maximal_positive_root_count_fixed_b_h(
    C::QQMatrix, M::ZZMatrix, L::QQMatrix,
    b_spec::Union{Vector{Int},Vector{QQFieldElem}}, 
    h::Union{Vector{Int},Vector{QQFieldElem}}; 
    TropB::Union{TropicalVariety,Nothing}=nothing, 
    TropL::Union{TropicalLinearSpace,Nothing}=nothing,
    verbose::Bool=false
)

    n = nrows(M) #number of variables
    m = ncols(M) #number of parameters
    s = rank(C) #rank
    d = n-s #corank

    @req nrows(L) == d "L must have the same number of rows as the corank of C"
    @req length(b_spec) == d "b_spec must have same length as the number of rows of L"
    @req length(h) == m "h must have same length as the number of columns of M"

    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:m)

    # Tropicalize the binomial part of the modified system
    if isnothing(TropB)
        binomials = vcat([y[i]-prod(x.^M[:,i]) for i=1:m], [z[1]-1])
        TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
        verbose && @info "Tropical binomial variety computed"
    end


    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(b))
    @req check_genericity_of_specialization(Lb, b_spec) "b_spec must be generic"
    Lb_spec = evaluate.(Lb, Ref(b_spec))

    if isnothing(TropL)
        linear_part_matrix = block_diagonal_matrix([Lb_spec, C])
        kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
        TropL = tropical_linear_space(kernel_matrix)
        verbose && @info "Tropical linear space computed"
    end

    pts, _ = tropical_stable_intersection_linear_binomial(TropL, TropB, 
        perturbation=vcat(zeros(Int, n+1), h), with_multiplicities=false)

    # Count how many of the tropical points that are positive
    Ilin = ideal(R, C*y) + ideal(R, Lb_spec*vcat(x,z))
    return count(Oscar.is_initial_positive(Ilin, nu, p) for p in pts)
end



@doc raw"""
    lower_bound_of_maximal_positive_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix; 
    num_b_attempts::Int=5, num_h_attempts_per_b::Int=10, verbose::Bool=false)

Computes a lower bound on the maximal positive root count of the augmented vertically parametrized 
system given by the coefficient matrix `C`, the exponent matrix `M`, and the affine form matrix `L`.

The function randomly samples `num_b_attempts` choices of the constant terms, and
for each such choice `num_h_attempts_per_b` shifts of the tropicalized binomial variety 
in the space of auxillary variables in the modification.


"""
function lower_bound_of_maximal_positive_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix; 
    num_b_attempts::Int=5, 
    num_h_attempts_per_b::Int=10, 
    show_progress::Bool=true,
    verbose::Bool=false
)
    n = nrows(M) #number of variables
    m = ncols(M) #number of parameters
    s = rank(C) #rank
    d = n-s #corank

    @req nrows(L) == d "L must have the same number of rows as the corank of C"

    # Tropicalize the binomial part of the modified system
    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:m)
    binomials = vcat([y[i]-prod(x.^M[:,i]) for i=1:m], [z[1]-1])
    TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
    verbose && @info "Tropical binomial variety computed"
   
    # Try different choices of b and h
    # Keep track of the maximal positive root count found and associated b and h values
    # Todo: Make this interruptible!
    best_count = 0
    best_b = nothing 
    best_h = nothing
    progress = ProgressMeter.Progress(num_b_attempts; 
        dt=0.4, 
        desc="Trying parameter values...", 
        barlen=30,
        output = stdout,
        enabled = show_progress
    );
    for b_attempt=1:num_b_attempts
        b_spec = L*rand(1:1000, n)
        Lb_spec = hcat(L, -matrix(L*rand(1:1000, n)))
    
        # Tropicalize the linear part of the modified system
        linear_part_matrix = block_diagonal_matrix([Lb_spec, C])
        kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
        TropL = tropical_linear_space(kernel_matrix)
        verbose && @info "Tropical linear space computed"
    
        # Compute the stable intersection for different h values
        new_count = nothing 
        h = nothing
        for h_attempt = 1:num_h_attempts_per_b
            generic_perturbation = false
            while !generic_perturbation
                try
                    h = rand(1:1000, m)
                    new_count = lower_bound_of_maximal_positive_root_count_fixed_b_h(
                        C, M, L, b_spec, h; TropB=TropB, TropL=TropL, verbose=verbose
                    )
                    generic_perturbation = true
                catch err
                    if isa(err, ErrorException) && err.msg == "random direction not generic"
                        continue
                    else
                        error(err)
                    end
                end
            end

            # Update the current best count
            if new_count > best_count
                best_count = new_count
                best_b = b_spec
                best_h = h
            end
        end

        # Update the progress bar
        ProgressMeter.update!(progress, b_attempt; 
            showvalues = [
                ("Number of b attempts", "$(b_attempt) ($(num_b_attempts))"), 
                ("Current maximal count", best_count)
            ]
        )
    end
    return best_count, best_b, best_h
end

@doc raw"""
    lower_bound_of_maximal_positive_steady_state_count(rn::ReactionSystem; kwargs...)

    Computes a lower bound on the maximal number of isolated positive steady states 
    that a mass action network `rn` can have.
"""
lower_bound_of_maximal_positive_steady_state_count(rn::ReactionSystem; kwargs...) = 
    lower_bound_of_maximal_positive_root_count(augmented_vertical_system(rn)...; kwargs...)
