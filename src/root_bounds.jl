
@doc raw"""
generic_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix; 
    b_spec = nothing, check_transversality::Bool=true, verbose::Bool=false)

    Compute the generic root count of an augmented vertically parametrized system given by 
    the coefficient matrix `C`, the exponent matrix `M`, and the affine form matrix `L`.

    # Example
    ```doctest
    julia> C = matrix(QQ, [1 -1 -1]);

    julia> M = matrix(ZZ, [1 0 2; 0 1 1]);
 
    julia> L = matrix(QQ, [1 1]);

    julia> generic_root_count(C, M, L, check_transversality=true)
    3

    julia> generic_root_count(C, M, L, check_transversality=false)
    3
    ```

"""
function generic_root_count(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M)); 
        b_spec = nothing, 
        check_transversality::Bool=true, 
        verbose::Bool=false)

    n = nrows(M) #number of variables
    m = ncols(M) #number of parameters
    s = rank(C) #rank
    d = n-s #corank

    @req nrows(L) == d "The augmentation matrix L must have the same number of rows as the corank of the coefficient matrix"

    # Coefficient matrix for the augmentation of the system
    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(b))

    # Pick a generic specialization of teh constant terms
    if isnothing(b_spec)
        is_generic = false
        while !is_generic
            b_spec = L*rand(1:1000, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
        end
    end
    @req check_genericity_of_specialization(Lb, b_spec) "Choice of constant terms to be generic"
    Lb_spec = evaluate.(Lb, Ref(b_spec))

    # Compute the generic root count as a mixed volume if the linear part gives a transversal matroid
    if check_transversality
        tp_nonlinear = transversal_presentation(C)
        tp_affine = transversal_presentation(Lb_spec)
        if tp_nonlinear != false && tp_affine != false
            verbose && @info "Transversal presentations found"
            nonlinear_supports = [Matrix{Int}(M[:,indices]) for indices in tp_nonlinear]
            affine_supports = [ hcat([i in 1:n ? standard_vector(i, n) : zeros(Int, n) for i in indices]...) for indices in tp_affine]
            supports = vcat(nonlinear_supports, affine_supports)
            return mixed_volume(supports)
        end
    end

    # Tropicalize the linear part of the modified system
    linear_part_matrix = block_diagonal_matrix([Lb_spec, C])
    kernel_matrix = transpose(kernel(linear_part_matrix, side=:right))
    TropL = tropical_linear_space(kernel_matrix)
    verbose && @info "Tropical linear space computed"

    # Tropicalize the binomial part of the modified system
    K, t = rational_function_field(QQ,"t")
    nu = tropical_semiring_map(K,t)
    R, x, z, y = polynomial_ring(K, "x"=>1:n, "z"=>1:1, "y"=>1:m)
    binomials = vcat([y[i]-prod(x.^M[:,i]) for i=1:m], [z[1]-1])
    TropB = Oscar.tropical_variety_binomial(ideal(R, binomials), nu)
    verbose && @info "Tropical binomial variety computed"
 
    # Compute the stable intersection
    generic_perturbation = false
    while !generic_perturbation
        try
            pts, mults = tropical_stable_intersection_linear_binomial(TropL, TropB)
            generic_perturbation = true
            return sum(mults)
        catch err
            if isa(err, ErrorException) && err.msg == "random direction not generic"
                continue
            else
                error(err)
            end
        end
    end
end



@doc raw"""
    steady_state_degree(rn::ReactionSystem; kwargs...)

Compute the steady state degree (in the complex torus) of the steady state system of a mass action network `rn`.

# Example
```jldoctest
julia> rn = @reaction_network begin
    k1, X1 --> X2
    k2, X2 --> X1
    k3, 2*X1 + X2 --> 3*X1
end;

julia> steady_state_degree(rn)
3
````

"""
steady_state_degree(rn::ReactionSystem; kwargs...) = 
    generic_root_count(augmented_vertical_system(rn)...; kwargs...)



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
```doctest
julia> C = matrix(QQ, [1 -1 -1]);

julia> M = matrix(ZZ, [1 0 2; 0 1 1]);

julia> L = matrix(QQ, [1 1]);

julia> julia> h = [37,97,18];

julia> b = [71];

julia> lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h)
3

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
    return count(is_initial_positive(Ilin, nu, p) for p in pts)
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

    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(b))

   
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

        # Pick a generic choice of b and specialize the augmentation
        is_generic = false
        b_spec = nothing #to make it accessible outside the while loop
        while !is_generic
            b_spec = L*rand(1:1000, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
        end
        Lb_spec = evaluate.(Lb, Ref(b_spec))

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

        # Update the progres bar (todo: figure out the showvalues bug)
        ProgressMeter.update!(progress, b_attempt; 
            #showvalues = [("Number of b attempts", "$(b_attempt) ($(num_b_attempts))"), ("Current maximal count", best_count)]
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


function toric_root_bound(A::ZZMatrix, L::QQMatrix;
    b_spec::Union{Nothing,Vector{Int},Vector{QQFieldElem}}=nothing,
    check_transversality::Bool=true,
    verbose::Bool=false
)

    n = ncols(A) # number of variables
    d = nrows(L) # number of affine equations
    @req rank(A) == d "System needs to be effectively square"

    # Add a column corresponding to the homogenization variable
    A_extended = hcat(A, zero_matrix(ZZ,nrows(A),1))

    # Pick a generic choice of constant terms
    B, b = rational_function_field(QQ, "b"=>1:d)
    Lb = hcat(B.(L), -matrix(b))

    # Pick a generic specialization of the constant terms
    if isnothing(b_spec)
        is_generic = false
        while !is_generic
            b_spec = L*rand(1:1000, n)
            is_generic = check_genericity_of_specialization(Lb, b_spec)
        end
    end
    @req check_genericity_of_specialization(Lb, b_spec) "Choice of constant terms to be generic"
    Lb_spec = evaluate.(Lb, Ref(b_spec))

    # Check for transversality
    # Warning: This currently gives incorrect values!
    if check_transversality
        tp = transversal_presentation(Lb_spec)
        if tp != false
            verbose && @info "Transversal presentation found"
            supports = [Matrix{Int}(A_extended[:,indices]) for indices in tp]
            degA = prod(diagonal(snf(A)))
            return mixed_volume(supports)/degA
        end
    end

    R, x, z = polynomial_ring(QQ, "x"=>1:n, "z"=>1:1)

    # Tropicalize the affine linear space
    TropL = tropical_linear_space(ideal(Lb_spec*vcat(x,z)))
    verbose && @info "Tropical linear space computed"

    # Tropicalize the toric variety
    I_toric = toric_ideal(R, transpose(A_extended))
    Trop_toric = Oscar.tropical_variety_binomial(I_toric,tropical_semiring_map(QQ))
    verbose && @info "Tropicalization of toric variety computed"
   
    # Stable intersection
    generic_perturbation = false
    while !generic_perturbation
        try
            pts, mults = tropical_stable_intersection_linear_binomial(TropL, Trop_toric)
            generic_perturbation = true
            return sum(mults)
        catch err
            if isa(err, ErrorException) && err.msg == "random direction not generic"
                continue
            else
                error(err)
            end
        end
    end
end