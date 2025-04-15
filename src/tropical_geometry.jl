
@doc raw""""
    tropical_stable_intersection_linear_binomial(TropL::TropicalLinearSpace,TropB::TropicalVariety)

Specialized stable intersection function for a tropical linear space 
and a linear space (encoded as a tropicalization of a binomial variety).

The output is a vector of stable intersection points and a vector with the multiplicities of the points.

"""
function tropical_stable_intersection_linear_binomial(TropL::TropicalLinearSpace, TropB::TropicalVariety; 
    perturbation::Union{Nothing, Vector{Int}} = nothing, with_multiplicities::Bool = true)

    bergmanRays, bergmanLineality = rays_modulo_lineality(TropL)
    bergmanRays = matrix(QQ, bergmanRays)
    bergmanLineality = matrix(QQ, bergmanLineality)

    minimalFaces, linearSpaceBasis = minimal_faces(TropB)
    linearSpaceBasis = matrix(QQ, linearSpaceBasis)

    @req length(minimalFaces) == 1 "Several minimal faces found in TropL"
    #perturbation = Vector(minimalFaces[1])
    #pick a random perturbation
    if isnothing(perturbation)
        perturbation = [-rand(-1000:1000) for i in 1:length(Vector(minimalFaces[1]))]
    end

    # compute the projection matrix onto the orthogonal complement of the euclidean linear space
    basisOfComplementTransposed = kernel(linearSpaceBasis, side=:right)
    basisOfComplement = transpose(basisOfComplementTransposed)
    projectionMatrix = basisOfComplementTransposed * inv(basisOfComplement * basisOfComplementTransposed) * basisOfComplement

    # project the rays of the Bergman fan
    projectedRays = bergmanRays * projectionMatrix
    projectedLineality = bergmanLineality * projectionMatrix

    #todo: add the bergmanLineality to the linearSpaceBasis

    # make it consistent whether projectionPerturbation and perturbation are rows/colums
    projectedPerturbation = matrix(QQ, [perturbation]) * projectionMatrix
    stableIntersectionPoints = Vector{QQFieldElem}[]
    stableIntersectionMults = Int[]

    indicesOfCones = ray_indices(maximal_polyhedra(TropL))
    nRaysPerCone = sum(indicesOfCones[1, :])
    for i in 1:nrows(indicesOfCones)
        # read off rays of the projected cone
        indicesOfCone = findall(indicesOfCones[i, :])
        projectedRaysOfCone = projectedRays[indicesOfCone, :]

        # test whether projected direction lies in projected cone
        # warning: be careful about the sign of the perturbation
        can_solve, solution = can_solve_with_solution(vcat(projectedRaysOfCone, projectedLineality),
            projectedPerturbation; side=:left)
        if can_solve
            firstZero = findfirst(isequal(0), solution)
            if (firstZero != nothing) && (firstZero[2] <= nRaysPerCone)
                # random direction lies on the boundary of the cone
                error("random direction not generic")
            end
            firstNegative = findfirst(a -> (a < 0), solution)
            if (firstNegative == nothing) || (firstNegative[2] > nRaysPerCone)
                # random direction lies in the interior of the cone,
                # compute intersection point and intersection multiplicity
                intersectionPoint = solution * vcat(bergmanRays[indicesOfCone, :], bergmanLineality)

                push!(stableIntersectionPoints, intersectionPoint[1, :])
                coneSpanBasis = vcat(bergmanRays[indicesOfCone, :], bergmanLineality)
                if with_multiplicities == true
                    push!(stableIntersectionMults, tropical_intersection_multiplicity(coneSpanBasis, linearSpaceBasis))
                end
            end
        end
    end
 
    return stableIntersectionPoints, stableIntersectionMults
end


function tropical_intersection_multiplicity(B1, B2)
    @assert ncols(B1) == ncols(B2) && nrows(B1) + nrows(B2) >= ncols(B1)

    # primitive scales every row by the lcm of the denominators, making the matrix integral
    # saturate computes a basis of the saturation of the sublattice spanned by the row vectors
    B1 = saturate(matrix(ZZ, Polymake.common.primitive(B1)))
    B2 = saturate(matrix(ZZ, Polymake.common.primitive(B2)))

    snfB12 = snf(vcat(B1, B2))
    return abs(prod([snfB12[i, i] for i in 1:ncols(snfB12)]))
end


@doc raw"""
    is_initial_positive(I::MPolyIdeal, nu::TropicalSemiringMap, w::AbstractVector)

Check whether a tropical point `w` is positive for the ideal `I` with respect to the tropical semiring map `nu`.
"""
function is_initial_positive(I::MPolyIdeal, nu::TropicalSemiringMap, w::AbstractVector)
    inI = initial(I, nu, w)
    G = Oscar.groebner_basis(inI; complete_reduction=true)
    @req all(isequal(2),length.(G)) "Gröbner basis must be binomial"
    # Check that the binomials have alternating signs
    return all(isequal(-1),[prod([sign(c) for c in Oscar.coefficients(g)]) for g in G])
end

