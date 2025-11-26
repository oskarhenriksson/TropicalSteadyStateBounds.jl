using TropicalSteadyStateBounds

for k = 1:3
    println("\n\n$(k)-site phosphorylation network:\n")
    C, M, L, A = multisite_phosphorylation_matrices(k)
    # t_grc = @elapsed grc = generic_root_count(C, M, L, check_transversality=false)
    # println("Generic root count: $grc (computed in $t_grc seconds)")
    t_toric = @elapsed grc_toric = toric_root_bound(A, L, check_transversality=false)
    println("Toric steady state bound: $grc_toric (computed in $t_toric seconds)")
    t_grc_mat = @elapsed grc_mat = generic_root_count(C, M, L, check_transversality=true)
    println("Generic root count with transversality: $grc_mat (computed in $t_grc_mat seconds)")
    t_toric_mat = @elapsed grc_toric_mat = toric_root_bound(A, L, check_transversality=true)
    println("Toric steady state bound with transversality: $grc_toric_mat (computed in $t_toric_mat seconds)")
end


# Old results before implementing the monomial reembedding
# The computation of the grc without transeversality only terminated for k=1

# 1-site phosphorylation network:

# Toric steady state bound: 3 (computed in 0.020498 seconds)
# Generic root count with transversality: 3 (computed in 0.012459917 seconds)
# Toric steady state bound with transversality: 3 (computed in 0.007808291 seconds)


# 2-site phosphorylation network:

# Toric steady state bound: 5 (computed in 0.035510417 seconds)
# Generic root count with transversality: 5 (computed in 0.016901917 seconds)
# Toric steady state bound with transversality: 5 (computed in 0.007522458 seconds)


# 3-site phosphorylation network:

# Toric steady state bound: 7 (computed in 0.22640075 seconds)
# Generic root count with transversality: 7 (computed in 0.591085084 seconds)
# Toric steady state bound with transversality: 7 (computed in 0.011513334 seconds)

