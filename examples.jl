using Catalyst
using TropicalSteadyStateBounds

# The running example
rn = @reaction_network begin
    k1, X1 --> X2
    k2, X2 --> X1
    k3, 2*X1 + X2 --> 3*X1
end;

# Compute the steady state degree directy from the network
@time sd = steady_state_degree(rn)

# The defining matrices for the steady state system
C, M, L = augmented_vertical_system(rn)

# Generic root count of the steady state system
@time generic_root_count(C, M, L)

# Without the transversality check
@time generic_root_count(C, M, L)

# Lower bound directly from the network
@time bound, b, h = lower_bound_of_maximal_positive_steady_state_count(rn)

@time bound, b, h = lower_bound_of_maximal_positive_root_count(C, M, L)

# Vertify the result for a given choice of b and h
h = [37,97,18]
b = [71]
lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h)


# Cell cycle
rn = Catalyst.@reaction_network begin
    k1, C + Mp --> C + M
    k2, Cp + M --> C + M
    k3, M + W --> Mp + W
    k4, M + W --> M + Wp
    k5, C --> Cp
    k6, Wp --> W
end;
@time steady_state_degree(rn)
@time lower_bound_of_maximal_positive_steady_state_count(rn)

# Verify the result for a given choice of b and h
C, M, L = augmented_vertical_system(rn)
h = [12, 86, 11, 27, 84, 98]
b =  [69, 42, 81]
lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h)

# HHK
rn = Catalyst.@reaction_network begin
    k1, HK00 --> HKp0
    k2, HKp0 -->  HK0p
    k3, HK0p --> HKpp  
    k4, HK0p  + Hpt --> HK00 + Hptp
    k5, HKpp  + Hpt --> HKp0 + Hptp
    k6, Hptp  --> Hpt
end;
@time steady_state_degree(rn)
@time lower_bound_of_maximal_positive_steady_state_count(rn)

# Verify the result for a given choice of b and h
C, M, L = augmented_vertical_system(rn)
h = [84, 46, 30, 13, 23, 68]
b =   [59, 34]
lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h)


# 1-site phosphorylation
rn = Catalyst.@reaction_network begin
  k1, S0 + E --> ES0
  k2, ES0  --> S0 + E
  k3, ES0  --> S1+E
  k4, S1 + F  --> FS1
  k5, FS1  --> S1 + F
  k6, FS1 --> S0 + F
end;

@time steady_state_degree(rn)
@time lower_bound_of_maximal_positive_steady_state_count(rn)

# Verify the result for a given choice of b and h
C, M, L = augmented_vertical_system(rn)
h =  [79, 26, 89, 92, 34, 83]
b =     [68, 52, 99]
lower_bound_of_maximal_positive_root_count_fixed_b_h(C, M, L, b, h)

#2-site phosphorylation
rn = Catalyst.@reaction_network begin
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

@time steady_state_degree(rn)





