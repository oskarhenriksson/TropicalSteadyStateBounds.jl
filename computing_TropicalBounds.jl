using Oscar
using Catalyst
include("TropicalBounds.jl")
include("tropical_geometry.jl")

#Running example
rn = @reaction_network begin
    k1, X1 --> X2
    k2, X2 --> X1
    k3, 2*X1 + X2 --> 3*X1
end 

C,M,L = get_augm_vert_param_system(rn)
@time generic_root_count(C,M,L)
p,h,b = maximum_positive_root_count_random_b(C,M,L,num_h_shifts= 8,num_random_b=5)
#After picking h and b randomly, we find
h = [37,97,18]
b = [-71]
maximum_positive_root_count(C,M,L,b,h)

#cell cycle
rn = @reaction_network begin
    k1, C + Mp --> C + M
    k2, Cp + M --> C + M
    k3, M + W --> Mp + W
    k4, M + W --> M + Wp
    k5, C --> Cp
    k6, Wp --> W
end
C,M,L = get_augm_vert_param_system(rn)
@time generic_root_count(C,M,L)
@time p,h,b = maximum_positive_root_count_random_b(C,M,L,num_h_shifts= 8,num_random_b=5)
#After picking h and b randomly, we find
h = [12, 86, 11, 27, 84, 98]
b =  [-69, -42, -81]
@time maximum_positive_root_count(C,M,L,b,h)

#HHK
rn = @reaction_network begin
    k1, HK00 --> HKp0
    k2, HKp0 -->  HK0p
    k3, HK0p --> HKpp  
    k4, HK0p  + Hpt --> HK00 + Hptp
    k5, HKpp  + Hpt --> HKp0 + Hptp
    k6, Hptp  --> Hpt
end
C,M,L = get_augm_vert_param_system(rn)
@time generic_root_count(C,M,L)
@time p,h,b = maximum_positive_root_count_random_b(C,M,L,num_h_shifts= 8,num_random_b=5)
#After picking h and b randomly, we find
h = [84, 46, 30, 13, 23, 68]
b =   [-59, -34]
@time maximum_positive_root_count(C,M,L,b,h)


#1-site phosphorylation
rn = @reaction_network begin
  k1, S0 + E --> ES0
  k2, ES0  --> S0 + E
  k3, ES0  --> S1+E
  k4, S1 + F  --> FS1
  k5, FS1  --> S1 + F
  k6, FS1 --> S0 + F
end 
C,M,L = get_augm_vert_param_system(rn)
@time generic_root_count(C,M,L)
@time p,h,b = maximum_positive_root_count_random_b(C,M,L,num_h_shifts= 8,num_random_b=5)
#After picking h and b randomly, we find
h =  [79, 26, 89, 92, 34, 83]
b =     [-68, -52, -99]
@time maximum_positive_root_count(C,M,L,b,h)

#2-site phosphorylation
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
C,M,L = get_augm_vert_param_system(rn)
@time generic_root_count(C,M,L)
@time p,h,b = maximum_positive_root_count_random_b(C,M,L,num_h_shifts= 10,num_random_b=2)
#After picking h and b randomly, we find
h = 
b =   
maximum_positive_root_count(C,M,L,b,h)








