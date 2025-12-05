using TropicalSteadyStateBounds

for k = 1:5
    println("\n\n$(k)-site phosphorylation network:\n")
    C, M, L, A = multisite_phosphorylation_matrices(k)
    if k<=2
        t_grc = @elapsed grc = generic_root_count(C, M, L, check_transversality=false)
        println("Generic root count: $grc (computed in $t_grc seconds)")
        t_pos = @elapsed positive_bound = lower_bound_of_maximal_positive_root_count(C, M, L; show_progress=false)[1]
        println("Lower bound of maximal positive steady state count: $positive_bound (computed in $t_pos seconds)")
    else
        println("Generic root count: skipped")
        println("Lower bound of maximal positive steady state count: skipped")
    end
    t_toric = @elapsed grc_toric = toric_root_bound(A, L, check_transversality=false)
    println("Toric steady state bound: $grc_toric (computed in $t_toric seconds)")
    t_toric_pos = @elapsed toric_positive_bound = toric_lower_bound_of_maximal_positive_root_count(A, L)[1]
    println("Toric lower bound of maximal positive steady state count: $toric_positive_bound (computed in $t_toric_pos seconds)")
    t_grc_mat = @elapsed grc_mat = generic_root_count(C, M, L, check_transversality=true)
    println("Generic root count with transversality: $grc_mat (computed in $t_grc_mat seconds)")
    t_toric_mat = @elapsed grc_toric_mat = toric_root_bound(A, L, check_transversality=true)
    println("Toric steady state bound with transversality: $grc_toric_mat (computed in $t_toric_mat seconds)")
end