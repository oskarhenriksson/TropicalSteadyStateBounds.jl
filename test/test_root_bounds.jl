using Test
using Random, Oscar, Catalyst

@testset verbose=true "Root bounds for networks" begin
@testset "Small example" begin

    C = matrix(QQ, [[1,-1,-1]])
    M = matrix(ZZ, [[1,0,2], [0,1,1]])
    L = matrix(QQ, [[1,1]])

    @test generic_root_count(C, M, L) == 3
    
    rn = @reaction_network begin
        k1, X1 --> X2
        k2, X2 --> X1
        k3, 2*X1 + X2 --> 3*X1
    end;

    @test steady_state_degree(rn, check_transversality=true) == 3
    @test steady_state_degree(rn, check_transversality=false) == 3

    h = [37,97,18]
    b = [71]
    @test lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h) == 3

    Random.seed!(1234)
    bound, _, _ = lower_bound_of_maximal_positive_steady_state_count(rn, num_b_attempts=5, num_h_attempts_per_b=3)
    @test bound == 3
    
end

@testset "Cell cycle" begin

    rn = @reaction_network begin
        k1, C + Mp --> C + M
        k2, Cp + M --> C + M
        k3, M + W --> Mp + W
        k4, M + W --> M + Wp
        k5, C --> Cp
        k6, Wp --> W
    end

    @test steady_state_degree(rn) == 2

    C, M, L = augmented_vertical_system(rn)
    h = [12, 86, 11, 27, 84, 98]
    b =  [69, 42, 81]
    @test lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h) == 2

end

@testset "HHK network" begin
    
    rn = @reaction_network begin
        k1, HK00 --> HKp0
        k2, HKp0 -->  HK0p
        k3, HK0p --> HKpp  
        k4, HK0p  + Hpt --> HK00 + Hptp
        k5, HKpp  + Hpt --> HKp0 + Hptp
        k6, Hptp  --> Hpt
    end

    @test steady_state_degree(rn) == 3

    C, M, L = augmented_vertical_system(rn)
    h = [84, 46, 30, 13, 23, 68]
    b =   [59, 34]
    @test lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h) == 3

end

@testset "1-site phosphorylation" begin

    rn = @reaction_network begin
        k1, S0 + E --> ES0
        k2, ES0  --> S0 + E
        k3, ES0  --> S1+E
        k4, S1 + F  --> FS1
        k5, FS1  --> S1 + F
        k6, FS1 --> S0 + F
    end

    @test steady_state_degree(rn) == 3

    C, M, L = augmented_vertical_system(rn)
    h =  [79, 26, 89, 92, 34, 83]
    b =     [68, 52, 99]
    @test lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h) == 1

end

@testset "2-site phosphorylation" begin

    rn = @reaction_network begin
    @parameters k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12
    @species E(t) F(t)  S0(t) S1(t) ES0(t) FS1(t) S2(t) ES1(t) FS2(t)
        k1, S0 + E --> ES0
        k2, ES0  --> S0 + E
        k3, ES0  --> S1+E
        k4, S1 + F  --> FS1
        k5, FS1  --> S1 + F
        k6, FS1 --> S0 + F
        k7, S1 + E --> ES1
        k8, ES1 --> S1 + E
        k9, ES1 --> S2 + E
        k10, S2 + F  -->FS2
        k11, FS2 --> S2 + F
        k12, FS2 --> S1 + F
    end 

    @test steady_state_degree(rn) == 5

end


@testset "Triangle network" begin


    rn = Catalyst.@reaction_network begin
        k1, 3*X1 + 2*X2 --> 6*X1
        k2, 3*X1 + 2*X2 --> 4*X2 
        k3, 4*X2 --> 3*X1 + 2*X2 
        k4, 6*X1 -->  4*X2
    end;

    C, M, L = augmented_vertical_system(rn)

    @test generic_root_count(C, M, L) == 6
    bound, _, _ = lower_bound_of_maximal_positive_root_count(C, M, L)
    @test bound == 1
    A = matrix(ZZ, [[3, 2]])
    @test toric_root_bound(A, L) == 3
    bound, _, _ = toric_lower_bound_of_maximal_positive_root_count(A, L)
    @test bound == 1

end

end